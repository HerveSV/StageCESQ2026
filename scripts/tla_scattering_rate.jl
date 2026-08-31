# ## Measure the photon scattering rate of a two-level atom for various illumination frequencies
using AtomTwin
using AtomTwin.Units
using Plots
using StatsBase
using Dates

# speed of light in vacuum, Boltzmann constant
#import PhysicalConstants.CODATA2022

# Fitting
using Optimization
using OptimizationOptimJL

# Convenience
using PyFormattedStrings


LOG_MESSAGES = true

HAVE_TWEEZER = true

# of Monte-Carlo shots for each detuning
N_shots = 500# 500#2000

s0_list = [10, 100]#[0.1, 10, 100]


# ## Parameters
const c_0 = 2.99792458e8
const k_B = 1.38064852e-23


# Falconi et al. PHYSICAL REVIEW LETTERS 135, 203402 (2025)
const T_init = 20μK     # Initial atom temperature

const waist_tweezer = 580nm

# We want to trap atoms at magic wavelength, where differential lightshift vanishes
# The Falconi paper provides the following value
const λ_magic_bluemot = 530nm
# Temperature equivalent of trap depth
const T_depth = 10mK#2.27mK           # 0


# Yb171 Blue MOT transition 1S0 -> 1P1
const Γ_bluemot = 2π * 29.13MHz
const ω_bluemot = 2π * 751.5274711THz

# Yb171 Yellow Clock transition 1S0 -> 3P0
const Γ_clock = 2π * 7.6e-3 
const ω_clock = 2π * 518.29583659086363THz
const λ_magic_clock = float(c_0) / 394.79947535THz



# Shorthand for the 
Γ = Γ_clock # Γ_bluemot
ω0 = ω_clock
λ_tweezer = 700nm#λ_magic_clock





# relation between Rabi frequency and saturation parameter for two-level atom
function tla_rabifreq(s0, Γ)
    return Γ * sqrt(s0/2)
end

function compute_FWHM(X, Y; max_idx=nothing)
    if isnothing(max_idx)
        y_max, max_idx = findmax(Y)
    else
        y_max = Y[max_idx]
    end
    y_half = y_max / 2.0

    # Interpolation function to find exact x where y = y_half between two indices
    interp_x(idx1, idx2) = X[idx1] + (y_half - Y[idx1]) * (X[idx2] - X[idx1]) / (Y[idx2] - Y[idx1])

    # Left side crossing (search backward from peak)
    left_idx2 = findlast(i -> Y[i] <= y_half, 1:max_idx)
    left_idx1 = left_idx2 + 1
    x_left = interp_x(left_idx1, left_idx2)

    # Right side crossing (search forward from peak)
    right_idx1 = max_idx + findfirst(i -> Y[i] <= y_half, (max_idx + 1):length(Y)) - 1
    right_idx2 = right_idx1 + 1
    x_right = interp_x(right_idx1, right_idx2)

    # FWHM calculation
    fwhm = x_right - x_left

    return fwhm

end


# ## System Definition
#
#
# Blue mot transition
#g, e = Level("1S0"), Level("1P1")
# Yellow clock transition
g, e = Level("1S0"), Level("3P0")

atom = Ytterbium171Atom(;
    levels = [g, e],
    x_init = [0.0, 0.0, 0.0],
    v_init = maxwellboltzmann(T=T_init)
)

if LOG_MESSAGES
    display(atom)
end


# ## Define a set of laser detunings
#
#
# Sharply peaked resonance, makes sense to sample more near resonance
sampling_method = "linear"
Δ_max = 7Γ                      #2π * 100e-3 # 100 kHz
Δ_min = 3e-4 * Γ                 #2π * 0.01e-3 #0.01MHz    # used for logarithmic sampling
γ = Γ * 2                        # width of dense sample area - for tangent sampling
N_half = 200                     # expect 2N_half + 1 points




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
elseif sampling_method == "test"
    detunings = 2π * 1e-3 .* range(82.5, 83.3, length=2N_half+1)#range(75.7, 100.0, length=100)
else
    Δ_step = Δ_max/N_points 
    detunings = range(-Δ_max, Δ_max, step=Δ_step)
end

if LOG_MESSAGES
    println(f"Detuning samples = {length(detunings)}, method = {sampling_method}")
    plt1 =Plots.scatter(detunings, detunings; title=f"samples = {length(detunings)}, method = {sampling_method}")
    display(plt1)
end


## Linear fit to extract scattering rate
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

beam = PlanarBeam(399e-9, 1.0, [0.0, 1.0, 0.0], [1, 0, 0])


# Define tunable parameters(can be changed at run-time)
Δ_param = Parameter(:delta, 0)
Ω_param = Parameter(:Omega, 1MHz)



if HAVE_TWEEZER    
    # Find power needed to achieve wanted trap
    U_depth = k_B * T_depth     # J
    # light shift per unit intensity for blue MOT transition
    U_I = AtomTwin._calc_light_shift(ω0, Γ, 2π*c_0/λ_tweezer)

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
else
    # Build full system without tweezer, just the two-level atom
    #system = System(atom)
end


# ## Drive two-level atom
#
#
initialize!(atom; beams=[tweezer])
coupling = add_coupling!(system, atom, g => e, Ω_param; active = false)
detuning = add_detuning!(system, atom, e, Δ_param; active = false)
lightshifts = add_lightshifts!(system, active = true)

# ## Build Sequence
#
# We want to measure emitted photons
# Also measure excited state population as sanity check
pd = PhotoDetectorSpec(name="click")
add_detector!(system, PhotoDetectorSpec(name="click"))
# need to register decay channel after detector is created
add_decay!(system, atom, e => g, Γ; clicks=pd)
add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))



# ## Run and analyse
# run a seperate simulation for each detuning

#pulse_duration =  3000 / (tla_rabifreq(1, Γ)/2π)
#dt = pulse_duration / 10_000

rates_list = []

if LOG_MESSAGES
    if HAVE_TWEEZER
        println(f"Tweezer wavelength: {λ_tweezer*1e9} nm")
    else
        println(f"No tweezer")
    end
end

# Create new directory to save data
compact_date = Dates.format(now(), "yyyymmdd_HHMMSS")
mkdir("./data/tla_scattering_rate/$compact_date")
for (j, s0) in enumerate(s0_list)
    Ω = tla_rabifreq(s0, Γ)
    scattering_rates = zeros(size(detunings))

    # adaptive time-steps
    pulse_duration =  100 / (Ω/2π)
    dt = 1e-2 / (Δ_max / 2π)#pulse_duration / 100_000

    if LOG_MESSAGES
        println(f"Saturation parameter s0: {s0}")
        println(f"dt = {dt:.2f}, pulse_duraction {pulse_duration:.2f}")
    end

    seq = Sequence(dt)
    @sequence seq begin
        Pulse([coupling, detuning], pulse_duration)
    end
    
    filepath = f"./data/tla_scattering_rate/{compact_date}/clock_s0={s0}_tweezer={λ_tweezer*1e9:.1f}nm_050826.csv"
    println(filepath)
    
    for (i, Δ_det) in enumerate(detunings)
        out = play(system, seq; initial_state=g, shots = N_shots, density_matrix = false, delta = Δ_det, Omega = Ω)

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

        if LOG_MESSAGES
            """
            plt0 = Plots.plot(
                tlist,
                mean_cumm_counts,
                xlabel    = "Time (s) [Trap on]",
                ylabel    = "Cummulative photon counts",
                color = :blue,
                linewidth = 1,
                alpha = 1,
                title = f"detuning = {detunings[i]/2π * 1e3} mHz, R/Γ = {rate/Γ}",

            )

            display(plt0)
            """
            # need to use (i-1)%10 as julia is 1-indexed
            if (i-1)%1 == 0
                println(f"Δ = {detunings[i]/2π * 1e3:.1f} mHz, R/Γ = {rate/Γ_clock:.2f}, lightshift = {system.nodes[3]._field._coeff}")
            end
        end
    end

    #fwhm = compute_FWHM(detunings, scattering_rates)

    plt_temp = Plots.plot(
        detunings./2π  .* 1e3,
        scattering_rates./Γ,
        xlabel = "Detuning (mHz)",
        ylabel = "R / Γ",
        label = f"s0 = {s0}",
        titlefont=font(12,"Arial"),
        title = f"s0={s0:.2f}"
    )

    display(plt_temp)

    #push!(fwhm_list, fwhm)
    push!(rates_list, scattering_rates)

    # ## Save data

    open(filepath, "w") do io
        println(io, "# Parameters: s0 = $s0, tweezer = $(λ_tweezer*1e9) nm, Rabi_freq = $Ω Hz, Gamma = $Γ Hz, N_shots = $N_shots, pulse_duration = $pulse_duration s, dt = $dt s")
       
        println(io, "Detuning [mHz],Scattering_Rate/Gamma")  # Your custom headers here
        for (det, rate) in zip(detunings./2π  .* 1e3, scattering_rates./Γ)
            println(io, "$det,$rate")
        end
    end
end

# Compute FWHM values
fwhm_list = []
for (j, s0) in enumerate(s0_list)
    fwhm = compute_FWHM(detunings, rates_list[j])

    push!(fwhm_list, fwhm)

end


# Prepare labels for plotting
label_list = []
for (j, s0) in enumerate(s0_list)
    label = f"{s0}, {fwhm_list[j]/2π * 1e-6:.1f} MHz, {Γ * sqrt(1+s0)/2π * 1e-6:.1f} MHz"
    push!(label_list, label)
end

print(size(rates_list))
alphas = [0, 0.5, 1]
colors = [:blue, :green, :red]
plt = Plots.plot(
    detunings./2π  .* 1e-6,
    rates_list./Γ,
    xlabel = "Detuning (MHz)",
    ylabel = "R / Γ",
    label = hcat(label_list...),
    ylim = [0, 0.5],
    #alpha = alphas,
    #color = colors
    titlefont=font(12,"Arial"),
    title = f"tweezer = {λ_tweezer*1e9:.1f}nm, Γ = {Γ/2π * 1e-6:.1f} MHz]"
)

Plots.hline!(plt, [0.5];
        alpha = 0.3,
        color = :red,
        linestyle = :dash,
        label = "",
    )

display(plt)

#savefig("./figs/notweezer_s0_scan_040826_2.png")
savefig(f"./data/tla_scattering_rate/{compact_date}/scatterrate_plot.pdf")
# we seem to get a linewidth which roughly follows Γ√(1+s0), which is the expected relation