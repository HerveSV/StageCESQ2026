using Dates

compact_str = Dates.format(now(), "yyyymmdd_HHMMSS")

mkdir("./data/tla_scattering_rate/$compact_str")