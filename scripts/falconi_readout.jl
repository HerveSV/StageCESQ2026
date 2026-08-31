# ## Number-resolved atom detection in tweezer trapped 171_Yb
# We aim to reproduce the results published in A. Muzi Falconi et al. Phys. Rev. Letters 135, 203402 (2025)

using AtomTwin
using AtomTwin.Units
using Plots


# ## Parameters

# Physical constants
const k_B = 1.380649e-23
const c = 29979245

# Tweezer
# PHYSICAL REVIEW LETTERS 135, 203402 (2025)
const waist = 580e-9            # Beam waist (m)
const λ_tweezer = 530e-9        # Tweezer wavelength (m)
const T_depth = 570e-6          # Trap depth, temperature equivalent U/k_B (K)

# 171_Yb
# M. Kroeze et al. 171Yb Reference Data (2026)
const g_I = -0.0005362
const g_L = 0.99999679
const g_S = 2.00231930436092

# PHYSICAL REVIEW LETTERS 135, 203402 (2025)
const temperature_init = 20e-6  # Intial temperature (K)

# no reference, randomly chosen
const Ω = 2π * 0.05e6           # Rabi frequency (rad/s)


# ## System Definition

function lande_J(J, L, S)
    return g_L * (J*(J+1) - S*(S+1) + L*(L+1)) / (2*J*(J+1)) + g_S * (J*(J+1) + S*(S+1) - L*(L+1)) / (2*J*(J+1))
end

function lande_F(F, J, I; g_J)
    #g_J = lande_J(J, L, S)
    return g_J * (F*(F+1) - I*(I+1) + J*(J+1)) / (2*F*(F+1)) + g_I * (F*(F + 1) + I*(I + 1) - J*(J + 1)) / ((2*F*(F+1)))
end


# F, J; label, g_F
# g_J values are experiments, taken from Kroeze et al. 171_Yb reference data
#ground = HyperfineManifold(1//2, 0; label="¹S₀", g_F=lande_F(1//2, 0, 1//2; g_J=0))
#excited_bluemot = HyperfineManifold(3//2, 0; label="¹P₁", g_F=lande_F(3//2, 0, 1//2; g_J=1.035))

# For now represent as two-level system without hyperfine structure
ground, excited_bluemot = Level("¹S₀"), Level("¹P₁")

#println(ground...)
yb = Ytterbium171Atom(; 
    levels=[ground, excited_bluemot],
    x_init=[0.0, 0.0, 0.0],
    v_init=maxwellboltzmann(T=temperature_init),
)


display(yb)


# Find the total power needed to achieve wanted trap depth 

U_depth = k_B * T_depth    # convert to Joules
# light shift for blue MOT transition
U_I = AtomTwin._calc_light_shift(2*π*751.527471e12, 29.13e6, 2*π*c/λ_tweezer)

# For fundamental Gaussian beam, peak intensity at waist is 
# related to total power via I0 = 2P/(πω₀²)
# Trap depth U0 = -U_I * I0, rearrange for P

P_tweezer = - 0.5 * π*(2*π*751.527471e12)^2 / (U_I)         # W


# Single-site tweezer array with specified geometry and powers
tweezer = GaussianBeam(
    λ    = λ_tweezer,
    w0   = waist,
    P   = P_tweezer         
)



# ## Test driving two-level atom

# Resonant planar beam driving the |g⟩ ↔ |e⟩ transition
beam = PlanarBeam(399e-9, 1.0, [1.0, 0.0, 0.0], [0, 1, 0])

# Build the full system and add a coherent coupling between |g⟩ and |e⟩
system   = System(yb, tweezer)

#coupling = add_coupling!(system, atom, g => e, Ω; beam = beam, active = false)
coupling = add_coupling!(system, yb, ground => excited_bluemot; active = false)

# ## Build Sequence
#
# We measure photon counts
# Also measure excited-state population for sanity check
add_detector!(system, PopulationDetectorSpec(yb, excited_bluemot; name = "P_e"))

PopulationDetectorSpec()