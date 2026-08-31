# ## Parameters

# Physical constants
c = 2.99792458e8

# PHYSICAL REVIEW LETTERS 135, 203402 (2025)
s0 = 40 # 1000           # s0 = I/I_Sat, for illumination beam
T_init = 20μK     # Initial atom temperature

waist_tweezer = 580nm
λ_tweezer = 730nm


# Yb171 Blue MOT transition 1S0 -> 1P1
const Γ_bluemot = 2π * 29.13MHz
const ω_bluemot = 2π * 751.5274711THz

# relation between Rabi frequency and saturation parameter for two-level atom
const Ω = Γ_bluemot * sqrt(s0/2)