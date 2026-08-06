using DataFrames, CSV
using Plots
using PyFormattedStrings

const Γ_bluemot = 2π * 29.13e6
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

# test effectiveness of FWHM finder
if false
    function lorentzian(x, x0, fwhm)
        z = 2(x-x0)/fwhm
        return 1/(1+z^2)
    end
w = 1.
x = range(-5., 5., 500)
y = lorentzian.(x, 0., w) # make lorenztian with FWHM=1

println(f"True value {w}")
println(f"Computed value {compute_FWHM(x, y)}")
end


s0_values = [0.1, 10.0, 100.0]

Δ_arrays = []
rate_arrays = []
FWHM_values = []
labels = []
labels_short = []

for (i, s0) in enumerate(s0_values)
    # 1. Read the CSV file into a DataFrame
    filepath = f"./data/tla_scattering_rate/s0={s0}_notweezer.csv"#f"./data/tla_scattering_rate/s0={s0}_tweezer=530.0nm.csv"
    df = CSV.read(filepath, DataFrame)

    # 2. Extract columns back into individual 1D arrays
    X = Vector(df[!,"Detuning [MHz]"])         # or df.X depending on your CSV headers
    Y = Vector(df[!,"Scattering_Rate/Gamma"])    # or df.Y

    FWHM = compute_FWHM(X, Y)

    push!(Δ_arrays, X)
    push!(rate_arrays, Y)
    push!(FWHM_values, FWHM)


    label = f"s0 = {s0}, FWHM = {FWHM:.2f} MHz, √2Ω = {tla_rabifreq(s0, Γ_bluemot)/2π * 1e-6:.2f} MHz"
    push!(labels, label)

    label_short = f"s0 = {s0}"
    push!(labels_short, label_short)
end

plt = Plots.plot(
    Δ_arrays[1], 
    rate_arrays,
    label = hcat(labels_short...),
    xlabel = "Detuning [MHz]",
    ylabel = "Scatter Rate / Γ",
    #title = "Tweezer wavelength: 730nm"
)
#Plots.scatter!(plt, [-FWHM/2, FWHM/2], [findmax(Y)[1]/2, findmax(Y)[1]/2])
display(plt)

#print(labels)

println(f"Γ = {Γ_bluemot/2π * 1e-6} MHz")
for l in labels
    println(l)
end


