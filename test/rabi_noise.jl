using AtomTwin
using StatsBase
using Plots


# frequency noise power spectral density (Hz2/Hz)
noise_model = LaserPhaseNoiseModel(
    bump_ampl = 2.0e5, # Hz2
    bump_center    = 1.0e6,
    bump_width     = 2e5,
    powerlaw_ampl  = 0.25e5,
    powerlaw_exp   = 0.0 # white noise
)

# is there any way I can plot this frequency PSD?
# yes, there is

# define a time interval
# collect() converts a range into an array
t = collect(range(0, 1e-3, step=1e-8)) # 1 ms with 1 ns steps
#print(t)

freqs, psd_values = AtomTwin.get_noise_spectrum(noise_model, t)

plot(freqs, psd_values, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="PSD (Hz²/Hz)", title="Laser Phase Noise Power Spectral Density")