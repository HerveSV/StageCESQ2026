import AtomTwin 

# ── Parameters ────────────────────────────────────────────────────────────────

Ω     = 2π * 1.0e6       # Rabi frequency (rad/s)
Γ     = 2π * 0.5e3       # Spontaneous decay rate |e⟩ → |g⟩ (rad/s)
T     = 1000 / (Ω/2π)    # Total time: 1000 Rabi periods (s)
dt    = T / 100_000      # 50 fixed steps per Rabi period
SHOTS = 100              # Monte Carlo trajectories
SAMPLES = 10



function build_atomtwin(; decay=false)
    g, e   = Level("g"), Level("e")
    atom   = Atom(; levels = [g, e])
    system = System(atom)
    decay && add_decay!(system, atom, e => g, Γ; active = true)
    add_detector!(system, PopulationDetectorSpec(atom, e; name = "P_e"))
    seq = Sequence(dt)
    @sequence seq begin Pulse(drive, T) end
    return system, seq, g
end

println("\n── AtomTwin ──────────────────────────────────────────────────────────")

let
    sys_u, seq_u, g = build_atomtwin(decay = false)
    sys,   seq,   _ = build_atomtwin(decay = true)

    job_se   = compile(sys_u, seq_u; initial_state = [g])
    job_me   = compile(sys,   seq;   density_matrix = true,  initial_state = [g])
    job_mcwf = compile(sys,   seq;   density_matrix = false, initial_state = [g])
    
    println("  Accuracy (max|err| last Rabi period, analytical reference):")
    report_acc("Schrödinger (unitary)",
        maxerr(play(job_se, sys_u).detectors["P_e"], P_se_ref))
    report_acc("master equation",
        maxerr(play(job_me, sys).detectors["P_e"], P_me_ref))

    P_mcwf = mean(play(job_mcwf, sys; shots = SHOTS).detectors["P_e"], dims = 2)
    report_acc("MCWF ($SHOTS shots avg)", maxerr(P_mcwf, P_me_ref))
    println()

    println("  Timing (minimum over samples):")
    report_time("Schrödinger (unitary)", @benchmark play($job_se,   $sys_u)            samples=SAMPLES evals=1)
    report_time("master equation",       @benchmark play($job_me,   $sys)              samples=SAMPLES evals=1)
    report_time("MCWF ($SHOTS shots, sequential)", @benchmark play($job_mcwf, $sys; shots=SHOTS, parallel_thresh=typemax(Int)) samples=SAMPLES evals=1)
    report_time("MCWF ($SHOTS shots, threaded)",   @benchmark play($job_mcwf, $sys; shots=SHOTS)                               samples=SAMPLES evals=1)
end