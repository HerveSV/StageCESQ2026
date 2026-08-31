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
const T_depth = 2.27mK           # 0


# Yb171 Blue MOT transition 1S0 -> 1P1
const Γ_bluemot = 2π * 29.13MHz
const ω_bluemot = 2π * 751.5274711THz

# Yb171 Yellow Clock transition 1S0 -> 3P0
const Γ_clock = 2π * 7.6e-3 
const ω_clock = 2π * 518.29583659086363THz
const λ_magic_clock = float(c_0) / 394.79947535THz


# Shorthand for the 
Γ = Γ_clock
ω0 = ω_clock
λ_tweezer = λ_magic_clock

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
    system = System(atom)
end


# ## Build Sequence

# We want to measure atom positions, to determine if the atom is in the trap
