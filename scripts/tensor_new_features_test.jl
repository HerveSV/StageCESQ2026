# Unit tests for tensor polarisability changes
using AtomTwin
using Plots
using CSV, DataFrames
#include("../src/atoms/polarizability.jl")
#include("../ext/AtomTwinPlotsExt.jl")

# ## 1. 1S0 and 3P0 curves
data_1S0 = CSV.read("/Users/hervesv/Documents/Stuff/Projects/AtomTwin.jl/examples/data/yb174/1S0.csv", DataFrame; delim=',', header=["Wavelength", "Lightshift"])

yb = Ytterbium171Atom(;)
models_0 = AtomTwin.getpolarizabilitymodels(yb)

print(models_0["1S0"])
typeof(models_0["1S0"])
curve = PolarizabilityCurve([models_0["1S0"], models_0["3P0"]])
fig1 = plot(curve)
scatter!(fig1, data_1S0[!, "Wavelength"], data_1S0[!, "Lightshift"], label="1S0 paper", alpha=0.5, markersize=3)


# ## 2. Defining non-Hyperfine manifold
# Expected: should throw warning "Tensor polarizability can only be assigned to Hyperfinelevel, '1S0' is of type Level; defaulting to α2 = 0.0"
g, e = Level("1S0"), Level("3P0")

atom = Ytterbium171Atom(;
    levels = [g, e]
)

tweezer = GaussianBeam(
        λ   = 500e-9,
        w0  = 1e-5,
        P   = 1.0
    )

initialize!(atom; beams=[tweezer])
system = System(atom, tweezer)

coupling = add_coupling!(system, atom, g => e, 1.0; active = false)
lightshifts = add_lightshifts!(system, atom)
seq = Sequence(1e-8)
@sequence seq begin
        Pulse([coupling], 1e-5)
end

out = play(system, seq)


# ## 3. Defining Hyperfine manifold