using AtomTwin
using Plots

# ── Parameters ────────────────────────────────────────────────────────────────

Ω               = 2π * 1.0e6        # Rabi frequency (rad/s)
Γ               = 2π * 0.5e3        # Spontaneous decay rate |e⟩ → |g⟩ (rad/s)
T               = 1000 / (Ω/2π)      # Total time: 100 Rabi periods (s)
dt              = T / 100_000       # 50 fixed steps per Rabi period

pulse_duration  = 30e-4            # Total pulse duration (s)
shots           = 1000              # Monte-Carlo shots

# 1. Create system
g, e   = Level("g"), Level("e")
atom   = Atom(; levels = [g, e])
system = System(atom)

coupling = add_coupling!(system, atom, g => e, Ω; active = false)
decay = add_decay!(system, atom, e => g, Γ; active=true)

# 2. Build system
add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))

seq = Sequence(dt)
@sequence seq begin
    Pulse(coupling, pulse_duration)
end

# 3. Run simulations

out_me = play(system, seq; initial_state = g, density_matrix = true) # master equation
out_qt = play(system, seq; initial_state = g, shots = 400) # quantum trajectories

# 4. Plot results

t = out_qt.times
t_us  = t .* 1e6   

Pe_me = out_me.detectors["P_e"]
Pe_qt = out_qt.detectors["P_e"]


plot(t_us, Pe_me)

