using AtomTwin
using Plots

## Default plot with inset


yb = Ytterbium171Atom(;)
models_0 = AtomTwin.getpolarizabilitymodels(yb)

curve = PolarizabilityCurve([models_0["1S0"], yb_models["3P0"]])
plot(curve)


# ## Need to redefine many structs

"""
    PolarizabilityModel_st

Empirical polarizability model for an atomic state, defined by a set of
discrete transitions and an optional offset.

# Fields
- `state::String`: Electronic state label (e.g. `"1S0"`, `"3P0"`).
- `transitions::Vector{NamedTuple}`: List of transitions; each entry has
  - `freq_THz::Float64`: Transition frequency in THz (linear). If negative, then state_e < state in energy.
  - `gamma_MHz::Float64`: Effective line width in MHz (linear) — see below.
  - `state_e::String`: Electronic configuration of excited state.
  - `Je::Rational`: Excited total angular momentum.
- `offset_Hz_per_Wm2::Float64`: Empirical offset in Hz/(W/m²).
- `reference::String`: Bibliographic reference for the data.


Note that transitions to states below the atomic state of interest (e.g. 1P1 -> 1S0, where we make a model for 1P1), 
`freq_THz` should be negative, in which case `Je` is interpreted as the ground total angular momentum and vice versa.
"""


struct PolarizabilityModel_st
    state::String
    Jg::Rational
    transitions::Vector{NamedTuple{(:freq_THz, :gamma_MHz, :state_e, :Je), Tuple{Float64, Float64, String, Rational}}}
    offset_Hz_per_Wm2::Float64
    reference::String
end

function PolarizabilityModel_st(state::String, Jg::Rational,
                             transitions::Vector;
                             offset_Hz_per_Wm2::Float64 = 0.0,
                             reference::String = "")
    norm = [_normalize_transition(t) for t in transitions]
    PolarizabilityModel(state, transit, offset_Hz_per_Wm2, reference)
end

const YB171_POLARIZABILITY_1S0_st = PolarizabilityModel_st(
    state = "1S0",
    Jg = 0,
    [
        (freq_THz = 539.386800, gamma_MHz = 0.183, state_e = "(6s6p) 3P1", Je = 1),     # (6s6p) 3P1
        (freq_THz = 751.526389, gamma_MHz = 29.127, state_e = "(6s6p) 1P1", Je = 1),    # (6s6p) 1P1
        (freq_THz = 865.111516, gamma_MHz = 11.052, state_e = "(7/2,5/2) J=1", Je = 1),    # (7/2,5/2) J=1
    ];
    offset_Hz_per_Wm2 = -0.8e-4,
    reference = "Phys. Rev. A 108, 053325 (2023)",
)





models_st = Dict(
    "1S0" => YB171_POLARIZABILITY_1S0,
    "3P0" => YB171_POLARIZABILITY_3P0,
    "3P1" => YB171_POLARIZABILITY_3P1
)
