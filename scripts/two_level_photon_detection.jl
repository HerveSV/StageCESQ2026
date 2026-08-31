# ## Number-resolved atom detection in tweezer trapped 171_Yb
# First step to reproducing the results published in A. Muzi Falconi et al. Phys. Rev. Letters 135, 203402 (2025)
# 
# We consider a tweezer trapped 2-level atom with resonant driving light
# We collect photons via photon detection

using AtomTwin
using AtomTwin.Units
using Plots
using StatsBase

# Fitting
using Optimization
using OptimizationOptimJL

# Convenicence
using PyFormattedStrings



# ## Parameters

# Physical constants
const c = 2.99792458e8
const k_B = 1.38064852e-23

# PHYSICAL REVIEW LETTERS 135, 203402 (2025)
const s0 = 40 # 1000           # s0 = I/I_Sat, for illumination beam
const T_init = 20μK     # Initial atom temperature

const waist_tweezer = 580nm
const λ_tweezer = 730nm
const T_depth = 5mK#2.27mK           # Temperature equivalent of trap depth


# Yb171 Blue MOT transition 1S0 -> 1P1
const Γ_bluemot = 2π * 29.13MHz
const ω_bluemot = 2π * 751.5274711THz

# relation between Rabi frequency and saturation parameter for two-level atom
const Ω = Γ_bluemot * sqrt(s0/2)
println(f"Rabi freq: {Ω/2π * 1e-6} MHz")
println(f"√2 Ω = {sqrt(2) * Ω/2π * 1e-6} MHz")
# ## System Definition

g, e = Level("1S0"), Level("3P0")

atom = Ytterbium171Atom(;
    levels = [g, e],
    x_init = [0.0, 0.0, 0.0],
    v_init = maxwellboltzmann(T=T_init)
)

display(atom)

# Find power needed to achieve wanted trap
U_depth = k_B * T_depth     # J
# light shift per unit intensity for blue MOT transition
U_I = AtomTwin._calc_light_shift(ω_bluemot, Γ_bluemot, 2π*c/λ_tweezer)

# For fundamental Gaussian beam, peak intensity I0 at waist is 
# related to total power P via I0 = 2P/(π w²)
# Trap depth U_depth = -U_I * I0, rearrange for P
#I0_tweezer = - U_depth / U_I
P_tweezer = - 0.5 * π * waist_tweezer^2 * U_depth/U_I

# Single-site tweezer
tweezer = GaussianBeam(
    λ   = λ_tweezer,
    w0  = waist_tweezer,
    P   = P_tweezer 
)



# ## Test driving two-level atom


# Resonant beam driving |g⟩ ↔ |e⟩ transition
# Propagates towards +y
# Horizontal polarisation along x
gaussian_illumination = false

if !gaussian_illumination
    # Plane wave if we want to address atom anywhere
    beam = PlanarBeam(399e-9, 1.0, [0.0, 1.0, 0.0], [1, 0, 0])
else
    # Gaussian beam if we want to address atoms only in trap
    beam = GaussianBeam(
        λ =  399nm,
        w0 = 4000nm, # 10 times the wavelength, collimated
        P = 1W # placeholder value
    )
end

#beam.E_field
#beam.unit_k
#beam.pol

# Build full system with coherent coupling between levels
system = System(atom, tweezer)

# Contructor for two level atom
coupling = add_coupling!(system, atom, g => e, Ω; active = false, beam = beam)
detuning = add_detuning!(system, atom, e, 0MHz; active = false)

initialize!(atom; beams=[tweezer])
lightshifts = add_lightshifts!(system, atom)

# ## Build Sequence
#
# We want to measure emitted photons
# Also measure excited state population as sanity check
# Measure atom positions to 

pd = PhotoDetectorSpec(name="click")
add_detector!(system, PhotoDetectorSpec(name="click"))
# need to register decay channel after detector is created
add_decay!(system, atom, e => g, Γ_bluemot; clicks=pd)

add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))

add_detector!(system, MotionDetectorSpec(atom; dims = [1, 2], name = "motion"))


pulse_duration = 10μs #1000 / (Ω/2π)
dt = 1e-10  #pulse_duration / 100_000

seq = Sequence(dt)
@sequence seq begin
    Pulse(coupling, pulse_duration)
end


# ## Run simulation
out = play(system, seq; initial_state=g, shots = 100, density_matrix = false)


# ## Analysis

tlist = out.times
counts = out.detectors["click"]
cumm_counts = cumsum(counts, dims=1)

# average over all shots
mean_cumm_counts = mean(cumm_counts, dims = 2)

# linear fit to extract scattering rate
linfunc(x; slope) = slope*x

# We must supply an objective function that will be minimized
# The u argument is a vector of parameters from the optimizer.
# data is a vector of static parameters passed through below.
function objective(u, data)
    # Get our fit parameters from u
    slope, = u
    # equivalent to:
    # slope = u[1]
    # intercept = u[2]

    # Get the x and y vectors from data
    x, y = data

    # Calculate the residuals between our model and the data
    residuals = linfunc.(x; slope) .- y

    # Return the sum of squares of the residuals to minimize
    return sum(residuals.^2)
end

# Define the initial parameter values for slope and intercept
u0 = [1.0]
# Pass through the data we want to fit
data = [tlist, mean_cumm_counts]

# Create an OptimizationProblem object to hold the function, initial
# values, and data.
prob = OptimizationProblem(objective, u0, data)

# Minimize the function. Optimization.jl uses the SciML common solver interface.
# Pass the problem you want to solve (optimization problem here) and a solver to use.
# NelderMead() is a derivative-free method for finding a function's local minimum.
sol = solve(prob, NelderMead())

# Exctract the best-fitting parameters
scatter_rate_data = sol.u[1]



scatter_rate_theo = Γ_bluemot * 0.5
cumm_counts_theo = tlist * scatter_rate_theo


println(f"Simulated scattering rate: {scatter_rate_data:.2f}")
println(f"Analytic scattering rate: {scatter_rate_theo:.2f}")

# ## Plotting

#=
# Plotting excited state population
plt = Plots.plot(
    tlist .* 1e6,
    out.detectors["P_e"],
    label     = "",
    xlabel    = "Time (μs)",
    ylabel    = "Population",
    title     = "Rabi oscillations with atomic motion",
    linewidth = 0.2,
    alpha     = 0.1,
    color     = :red,
    legend    = :none,
)


Plots.plot!(
    plt,
    tlist .* 1e6,
    mean(out.detectors["P_e"], dims = 2),
    color     = :black,
    linewidth = 3,
)

plt
=#


plt2 = Plots.plot(
    tlist .* 1e6,
    mean_cumm_counts,
    xlabel    = "Time (μs)",
    ylabel    = "Cummulative photon counts",
    color = :blue,
    linewidth = 1,
    alpha = 1,
    label = "Data",
)

Plots.plot!(
    plt2, 
    tlist .* 1e6,
    cumm_counts_theo,
    color     = :red,
    linewidth = 1,
    alpha = 0.8,
    label = "Theory",
    title = "730nm Lightshift on"
)

display(plt2)

#savefig("figs/759_nolightshift.pdf")

atom_motion = out.detectors["motion"]

# atom_motion has dims (time, channels, shots) e.g. 1000 x 2 x 100
# stack all shots sequentially to get (time*shots) x channels -> 100000 x 2
atom_motion = out.detectors["motion"]
nt = size(atom_motion, 1)
nshots = size(atom_motion, 3)
nch = size(atom_motion, 2)
atom_motion_stacked = reshape(permutedims(atom_motion, (1, 3, 2)), nt * nshots, nch)

plt3 = plot(
    collect(1:size(atom_motion_stacked, 1)),
    atom_motion_stacked[:, 1],
    xlabel = "Sample index",
    ylabel = "Motion channel 1",
    label = "All shots stacked"
)

# XY position heatmap from the motion detector
xpos = atom_motion_stacked[:, 1]
ypos = atom_motion_stacked[:, 2]

plt_xy = histogram2d(
    xpos,
    ypos,
    nbins = (80, 80),
    xlabel = "x position",
    ylabel = "y position",
    title = "Atom position XY heatmap",
    colorbar = true,
    aspect_ratio = 1,
    legend = false,
    c = :viridis,
)

display(plt_xy)


