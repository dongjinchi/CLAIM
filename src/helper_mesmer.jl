using CSVFiles, DataFrames

function load_region_lt(esm)
    path = "data/mesmer_data/emulations/local/region_T/"
    region_lt_coef  = DataFrame(load(joinpath(@__DIR__, "../"*path*esm*"_regional_T_coef.csv")))
    region_lt_inter = DataFrame(load(joinpath(@__DIR__, "../"*path*esm*"_regional_T_inter.csv")))
    return region_lt_coef, region_lt_inter
end

function load_region_lv(esm,emu)
    path = "data/mesmer_data/emulations/local/region_variability/"
    region_lv  = DataFrame(load(joinpath(@__DIR__, "../"*path*esm*"_regional_variability.csv")))
    region_lv = filter(row -> row.emu == emu, region_lv)
    select!(region_lv, Not(:emu))
    return region_lv
end
