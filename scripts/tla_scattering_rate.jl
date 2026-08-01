# ## Measure the photon scattering rate of a two-level atom for various illumination frequencies
using AtomTwin
using AtomTwin.Units
using Plots
using StatsBase
using DataFrames, CSV 

# Fitting
using Optimization
using OptimizationOptimJL

# Convenience
using PyFormattedStrings

# ## Parameters

# Physical constants
const k_B = 1.380649e-23
const c = 29979245

# Falconi et al. PHYSICAL REVIEW LETTERS 135, 203402 (2025)
const T_init = 20μK     # Initial atom temperature

const waist_tweezer = 580nm

# We want to trap atoms at magic wavelength, where differential lightshift vanishes
# The Falconi paper provides the following value
const λ_magic_bluemot = 530nm
# Temperature equivalent of trap depth
const T_depth = 2.27mK           


# Yb171 Blue MOT transition 1S0 -> 1P1
const Γ_bluemot = 2π * 29.13MHz
const ω_bluemot = 2π * 751.5274711THz

# Yb171 Yellow Clock transition 1S0 -> 3P0
const Γ_clock = 2π * 7.6e-3 
const ω_clock = 2π * 518.29583659086363THz


# Shorthand for the 
const Γ = Γ_bluemot
const ω0 = ω_bluemot
const λ_tweezer = λ_magic_bluemot


const pulse_duration = 0.1μs #1000 / (Ω/2π)
const dt = 1e-10#pulse_duration / 100_000

# relation between Rabi frequency and saturation parameter for two-level atom
function tla_rabifreq(s0, Γ)
    return Γ * sqrt(s0/2)
end


# ## System Definition
#
#

# Blue mot transition
g, e = Level("1S0"), Level("1P1")
# Yellow clock transition
#g, e = Level("1S0"), Level("3P0")

atom = Ytterbium171Atom(;
    levels = [g, e],
    x_init = [0.0, 0.0, 0.0],
    v_init = maxwellboltzmann(T=T_init)
)

display(atom)


#print(Γ/2π * 1e-6)
# ## Define a set of laser detunings
# Sharply peaked resonance, makes sense to sample more near resonance
sampling_method = "tan"
Δ_max = 2π * 200MHz
Δ_min = 2π * 0.01MHz    # used for logarithmic sampling
γ = Γ_bluemot * 2 # width of dense sample area - for 
N_half = 100 # expect 2N_half + 1 points


s0_list = [0.1, 10, 100]
#λ_list = [λ_tweezer-200nm, λ_tweezer, λ_tweezer+200nm] # around magic wavelength

if sampling_method == "log"
    # Generate points from Δ_min to Δ_max in log space
    pos_side = 10 .^ range(log10(Δ_min), log10(Δ_max), length=N_half)
    # Mirror to negative side and include exact center
    detunings = sort(unique(vcat(.- pos_side, 0., .+ pos_side)))
elseif sampling_method == "tan"
    Θ_max = atan(Δ_max / γ)
    Θ_grid = range(-Θ_max, Θ_max, length=2*N_half+1)
    detunings = γ .* tan.(Θ_grid)
elseif sampling_method == "linear"
    Δ_step = Δ_max/N_half
    detunings = range(-Δ_max, Δ_max, step=Δ_step)
else
    println("Default to linear sampling for detuning")
    Δ_step = Δ_max/N_points 
    detunings = range(-Δ_max, Δ_max, step=Δ_step)
end

if false
    Plots.scatter(detunings, detunings)
end

scattering_rates = zeros(size(detunings))

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

size(detunings)[1]

beam = PlanarBeam(399e-9, 1.0, [0.0, 1.0, 0.0], [1, 0, 0])


# Define tunable parameters(can be changed at run-time)
Δ_param = Parameter(:delta, 0)
Ω_param = Parameter(:Omega, 1MHz)


# Find power needed to achieve wanted trap
U_depth = k_B * T_depth     # J
# light shift per unit intensity for blue MOT transition
U_I = AtomTwin._calc_light_shift(ω0, Γ, 2π*c/λ_tweezer)

# For fundamental Gaussian beam, peak intensity I0 at waist is 
# related to total power P via I0 = 2P/(π w²)
# Trap depth U_depth = -U_I * I0, rearrange for P
P_tweezer = - 0.5 * π * waist_tweezer^2 * U_depth/U_I

# Single-site tweezer, tunable wavelength
tweezer = GaussianBeam(
    λ   = λ_tweezer,
    w0  = waist_tweezer,
    P   = P_tweezer 
)
# Build full system with coherent coupling between levels
system = System(atom, tweezer)



# ## Drive two-level atom
coupling = add_coupling!(system, atom, g => e, Ω_param; active = false, beam = beam)
detuning = add_detuning!(system, atom, e, Δ_param; active = false)

# ## Build Sequence
#
# We want to measure emitted photons
# Also measure excited state population as sanity check
pd = PhotoDetectorSpec(name="click")
add_detector!(system, PhotoDetectorSpec(name="click"))
# need to register decay channel after detector is created
add_decay!(system, atom, e => g, Γ; clicks=pd)
add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))

seq = Sequence(dt)
@sequence seq begin
    Pulse([coupling, detuning], pulse_duration)
end



# ## Run and analyse
# run a seperate simulation for each detuning

println(f"Tweezer wavelength: {λ_tweezer*1e9} nm")
for (j, s0) in enumerate(s0_list)
    println(f"Saturation parameter s0: {s0}")
    println(f"./data/tla_scattering_rate/s0={s0}_tweezer={λ_tweezer*1e9}nm.csv")
    Ω = tla_rabifreq(s0, Γ)
    for (i, Δ_det) in enumerate(detunings)
        out = play(system, seq; initial_state=g, shots = 5000, density_matrix = false, delta = Δ_det, Omega = Ω)

        tlist = out.times
        counts = out.detectors["click"]
        cumm_counts = cumsum(counts, dims=1)

        # average over all shots
        mean_cumm_counts = mean(cumm_counts, dims = 2)

        # Linear fit to extract scattering rate
        
        u0 = [1.0]  # Define the initial parameter values for slope
        data = [tlist, mean_cumm_counts]    # Pass through the data we want to fit
        # Create an OptimizationProblem object to hold the function, initial
        # values, and data.
        prob = OptimizationProblem(objective, u0, data)
        # Minimize the function
        sol = solve(prob, NelderMead())
        rate = sol.u[1]

        scattering_rates[i] = rate
        if i%10 == 0
            println(f"Δ = {detunings[i]/2π * 1e-6:.1f} MHz, R/Γ = {rate/Γ_bluemot:.2f}")
        end
    end

    # ## Plotting

    Plots.plot(
        detunings/2π  * 1e-6,
        scattering_rates/Γ,
        xlabel = "Detuning (MHz)",
        ylabel = "R / Γ",
        label = "Simulated values",
        ylim = [0, 0.5]
    )


    # ## Save data

    #CSV.write(f"../data/tla_scattering_rate/s={s}.csv", (X = X, Y = Y))

    open(f"./data/tla_scattering_rate/s0={s0}_tweezer={λ_tweezer*1e9}nm.csv", "w") do io
        println(io, "Detuning [MHz],Scattering_Rate/Gamma")  # Your custom headers here
        for (det, rate) in zip(detunings/2π  * 1e-6, scattering_rates/Γ_bluemot)
            println(io, "$det,$rate")
        end
    end
end

Plots.hline([0.5];
        alpha = 0.3,
        color = :red,
        linestyle = :dash,
        label = "",
    )