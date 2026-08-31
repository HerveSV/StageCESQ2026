
using WignerSymbols
using Plots
using Dates
using PyFormattedStrings
using CSV
using DataFrames

# Computes hyperfine reduced dipole matrix in the LS coupling approximation
function reduced_hyperfine_matrix(numbers, numbers_prime)
    l1, l2, L, S, J = numbers["l1"], numbers["l2"], numbers["L"], numbers["S"], numbers["J"]
    l1_p, l2_p, L_p, S_p, J_p = numbers_prime["l1"], numbers_prime["l2"], numbers_prime["L"], numbers_prime["S"], numbers_prime["J"]

    if S != S_p
        return 0
    elseif (l1 == l2) && (l1_p == l2_p)
        return 0
    end
end

df_even = CSV.read("./yb_even.csv", DataFrame; delim=';')
df_odd = CSV.read("./yb_odd.csv", DataFrame; delim=';')

# 3D1
top_energy = 24489.102
#df_even[2, "Energy[cm-1]"]

# we want to find all lower lying odd states

mask = zeros(size(df_odd)[1])
i = 1
for data in eachrow(df_odd)
    if data["Energy[cm-1]"] < top_energy
        mask[i] = 1
    end
    i += 1
end


mask = df_odd."Energy[cm-1]" .< top_energy

lower_levels = df_odd[mask, :]

println(lower_levels)

# Consider only the LS coupled levels

Levels_pp = ["3P0", "3P1", "3P2", "(7/2, 3/2)2"]
J_pp = [0, 1, 2]
L_pp = [1, 1, 1]
energies = [17288.4, 17992.0, 19710.4]

# Top level 3D1
J_p = 1
L_p = 2

# Desire bottom level 3P0
J = J_pp[1]
L = L_pp[1]
S = 1

# we are interest in the branching ratio 3D1 -> 3P0
numerator = (top_energy - energies[1])^3 * wigner6j(L, L_p, 1, J_p, J, S)^2

denom = (top_energy .- energies).^3 .* wigner6j.(L, L_p, 1, J_p, J_pp, S).^2
denom = sum(denom)

cumsum(denom) .- sum(denom)
sum(denom)

beta = numerator / denom
beta = numerator / denom

4.922610110144811252752939860026041666666666666666666666666666666666666666666627e+10 - 4.922610110144811252752939860026041666666666666666666666666666666666666666666627e+10

# Lets consider only the LS coupled ones


# Filter rows where Designation equals "3D"
#df_filtered = filter(row -> row["Level energy"] == 0, df_even)
#df_filtered = filter(row -> row["Configuration"] == "4f14(1S)6s6d", df_even)


function relevant_transitions(df_transitions, J_init, energy)

    mask1 = df_transitions."Energy[cm-1]" .< energy
    mask2 = (df_transitions."J" .== J_init) .|| (df_transitions."J" .== J_init+1) .|| (df_transitions."J" .== J-1)

    mask_tot = mask1 .& mask2

    levels = df_transitions[mask_tot, :]
    return levels
end


# THIS WORKS!!!
# What about (5d6s)3D2
energy_3D2 = 24751.948
J_3D2 = 2

levels2 = relevant_transitions(df_odd, J_3D2, energy_3D2)


Levels_pp = ["3P0", "3P1", "3P2", "(7/2,  3/2)2"]
J_pp = [0, 1, 2]#, 1]
L_pp = [1, 1, 1]#, 2]
energies = [17288.4, 17992.0, 19710.4]#, 23188.5]

# Top level 3D2
J_p = 2
L_p = 2

# Desire bottom level 3P1
J = 1
L = 1
S = 1

# we are interest in the branching ratio 3D1 -> 3P0
numerator = (energy_3D2 - energies[2])^3 * (2J + 1) *  wigner6j(L, L_p, 1, J_p, J, S)^2

denom = (energy_3D2 .- energies).^3 .* (2J_pp .+ 1) .* wigner6j.(L, L_p, 1, J_p, J_pp, S).^2
denom = sum(denom)


beta = numerator / denom


# This works too!!
# What about (6s5d)3D1
energy_3D1 = 24489.102
J_3D1 = 1
L_3D1 = 2
levels_3 = relevant_transitions(df_odd, J_3D1, energy_3D1)

Levels_pp = ["3P0", "3P1", "3P2"]
J_pp = [0, 1, 2]#, 1]
L_pp = [1, 1, 1]#, 2]
energies = [17288.4, 17992.0, 19710.4]#, 23188.5]

# Top level 3D1
J_p = 1
L_p = 2

# Desire bottom level 3P0
J = 0
L = 1
S = 1

numerator = (energy_3D1 - energies[1])^3 * (2J + 1) *  wigner6j(L, L_p, 1, J_p, J, S)^2

denom = (energy_3D1 .- energies).^3 .* (2J_pp .+ 1) .* wigner6j.(L, L_p, 1, J_p, J_pp, S).^2
denom = sum(denom)

beta = numerator / denom


# What about (6s7s)3S1
energy_3S1 = 32694.692 # cm-1
J_3S1 = 1

levels_1 = relevant_transitions(df_odd, J_3S1, energy_3S1)
print(levels_1)

# ignore none-jj-coupled states
Levels_pp = ["3P0", "3P1", "3P2", "1P1"]
J_pp = [0, 1, 2, 1]
L_pp = [1, 1, 1, 1]
energies = [17288.4, 17992.0, 19710.4, 25068.2]#, 23188.5]

# Top level 3S1
J_p = 1
L_p = 0

# Desire bottom level 3P1
J = 1
L = 1
S = 1

# we are interest in the branching ratio 3S1 -> 3P1
numerator = (energy_3S1 - energies[2])^3 * (2J + 1) *  wigner6j(L, L_p, 1, J_p, J, S)^2

denom = (energy_3S1 .- energies).^3 .* (2J_pp .+ 1) .* wigner6j.(L, L_p, 1, J_p, J_pp, S).^2
denom = sum(denom)

beta = numerator / denom


# we are interest in the branching ratio 3S1 -> 3P0
numerator = (energy_3S1 - energies[1])^3 * (2J + 1) *  wigner6j(L, L_p, 1, J_p, J, S)^2

denom = (energy_3S1 .- energies).^3 .* (2J_pp .+ 1) .* wigner6j.(L, L_p, 1, J_p, J_pp, S).^2
denom = sum(denom)

beta = numerator / denom