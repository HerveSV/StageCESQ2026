using AtomTwin
using Plots # Or Makie for visualizing the output

# Step 1: define the system
g, e = Level(; label = "g"), Level(; label = "e")
atom = Atom(; levels = [g, e])
system = System(atom)

# Coupling with Rabi frequency 1 MHz
add_coupling!(system, atom, g => e, 2pi * 1e6) # Omega/2pi = 1 MHz
# Add spontaneous emission
Γ = 2π * 0.3e6
add_decay!(system, atom, e => g, Γ)

add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))


# Step 2: build a sequence
seq = Sequence(1e-9) # fixed time step: 1 ns
@sequence seq begin
    Wait(5e-6) # evolve for 5 us
end
# Step 3: run and inspect
# if density_matrix=false while dissipative channels are present, will do monte-carle quantum jumps
out = play(system, seq; initial_state = g, density_matrxi=true)
# out.detectors["P_e"] is a Vector{Float64} of P_e(t)
# out.times is the corresponding time axis

# 5. Extract and plot results
times = out.times
pe = out.detectors["P_e"]

plot(times * 1e6, pe, 
     xlabel = "Time (μs)", 
     ylabel = "Excited State Population", 
     title = "2-Level Atom Rabi Oscillations",
     lw = 2, legend = false)