# ============================================================================
# 【helper_mesmer.jl】MESMER 网格参数加载、局地变率模拟与国家聚合。
#
# 设计背景（open_issues.md P1-3）
# ------------------------------------------------------------------
# 旧链路把 MESMER 的 5,121 个陆地格点先平均到 2,388 个次国家 RegionID，
# 再按面积加权聚合到 185 国，两层平均都是等网格/面积口径，与回归端
# lctmp_bkly_pw 的人口加权口径不一致；RegionID 映射还静默丢国家，导致
# 6 国 β=γ=0（前向模型里永不升温）、26 国 region_lv 恒零。
#
# 现在改为：mesmer_data/export_grid_params.py 只做训练并持久化网格级参数，
# 本文件在 Julia 端完成
#   1. 网格 β/γ  --W-->  国家 ltcoef/ltinter；
#   2. 网格局地变率模拟（AR(1) + 空间相关新息）--W--> 国家 region_lv。
# W 是 185×5940 的人口权重矩阵，由 03Claim/scripts/build_mesmer_country_weights.py
# 生成，与回归端共享同一份 GlobPOP、同一份 world_countries.shp 和同一张
# TEMPERATURE_PROXY 代理表。
#
# 为什么"先模拟后聚合"这个顺序不能颠倒
# ------------------------------------------------------------------
# MESMER 的局地变率是 AR(1) + 空间相关新息（OLS_AR1_sci）。若先把参数聚合
# 到国家再模拟：(a) 国家间的天气相关性会丢失，除非另行投影新息协方差；
# (b) 各格点 AR(1) 系数不同，其加权和是 ARMA 而非 AR(1)，自相关口径会变。
# 因此模拟必须在网格上完成，聚合放最后一步。
#
# 与 Python 端的可复现性
# ------------------------------------------------------------------
# numpy 的 multivariate_normal 默认用 SVD 分解协方差并使用其自身 RNG，
# Julia 端用 Cholesky + Julia RNG，因此**无法逐位复现** Python 的实现；
# 旧的 300 个 emu 文件不可能被重现。等价性只能在分布层面检验（格点方差、
# lag-1 自相关、空间相关结构、聚合后国家级方差与 ACF）。复现性改由本文件
# 的显式 seed 保证，seed 会写进输出文件名与 manifest。
#
# 旧的 load_region_lt / load_region_lv 读取 mesmer_results/ 下按 RegionID
# 聚合的 302 个存档文件，现保留在 §存档 一节，仅用于与新链路做新旧对比，
# 不参与模型运行。
# ============================================================================

using CSVFiles, DataFrames, DelimitedFiles, Random, SHA, SparseArrays, Statistics

# MESMER 内部的预热长度，见 mesmer/create_emulations/create_emus_gv.py:153
# 与 create_emus_lv.py:204。必须与上游一致，否则起始段方差偏低。
const MESMER_GV_BUFFER = 50
const MESMER_LV_BUFFER = 20

# 旧 aggregate_mesmer_variability 把 year <= 2024 的变率强制置零：历史段用的是
# 观测天气而非 emulator 抽样。保留该口径，否则历史损害路径会无声改变。
const MESMER_VARIABILITY_START_YEAR = 2025

# ----------------------------------------------------------------------------
# 一、网格参数加载
# ----------------------------------------------------------------------------

# 【_read_numeric_csv】读取带表头的纯数值 CSV。3.6M 行的 Cholesky 三元组用
# _read_csv_dicts 会非常慢，故这里走 DelimitedFiles。
function _read_numeric_csv(path::AbstractString)
    isfile(path) || error("missing MESMER grid parameter file: $path")
    data, header = readdlm(path, ',', Float64, '\n'; header=true)
    cols = Dict(strip(String(h)) => i for (i, h) in enumerate(vec(header)))
    return data, cols
end

function _grid_column(data, cols, name::AbstractString)
    haskey(cols, name) || error("column '$name' not found; have $(sort(collect(keys(cols))))")
    return data[:, cols[name]]
end

function _check_contiguous_land_idx(idx, file)
    idx == collect(0:length(idx)-1) ||
        error("$file: land_idx must be contiguous 0..n-1 in row order")
end

# 【load_grid_lt】网格级局地趋势参数 β（coef_gttas）与 γ（intercept）。
function load_grid_lt(unshared_dir::AbstractString)
    data, cols = _read_numeric_csv(joinpath(unshared_dir, "climate_grid_lt.csv"))
    _check_contiguous_land_idx(Int.(_grid_column(data, cols, "land_idx")), "climate_grid_lt.csv")
    return (beta=_grid_column(data, cols, "coef_gttas"),
            gamma=_grid_column(data, cols, "intercept"))
end

# 【load_grid_lv】网格级局地变率参数：全球变率载荷 coef_gvtas 与 AR(1) 三件套。
function load_grid_lv(unshared_dir::AbstractString)
    data, cols = _read_numeric_csv(joinpath(unshared_dir, "climate_grid_lv.csv"))
    _check_contiguous_land_idx(Int.(_grid_column(data, cols, "land_idx")), "climate_grid_lv.csv")
    return (coef_gvtas=_grid_column(data, cols, "coef_gvtas"),
            ar1_int=_grid_column(data, cols, "AR1_int"),
            ar1_coef=_grid_column(data, cols, "AR1_coef"),
            ar1_std=_grid_column(data, cols, "AR1_std_innovs"))
end

# 【load_grid_gv】全球变率 AR(p) 参数（当前标定为 AR(2)）。
function load_grid_gv(unshared_dir::AbstractString)
    data, cols = _read_numeric_csv(joinpath(unshared_dir, "climate_grid_gv.csv"))
    order = Int(data[1, cols["ar_order"]])
    coefs = [data[1, cols["ar_coef$(i)"]] for i in 1:order]
    return (order=order, ar_int=data[1, cols["ar_int"]],
            ar_std=data[1, cols["ar_std_innovs"]], coefs=coefs)
end

# 【load_innov_cholesky】新息协方差的下三角 Cholesky 因子（稀疏三元组）。
# 文件中的 row/col 为 0-based，转为 Julia 的 1-based。
function load_innov_cholesky(unshared_dir::AbstractString, n_grid::Int)
    data, cols = _read_numeric_csv(joinpath(unshared_dir, "climate_grid_innov_chol.csv"))
    rows = Int.(_grid_column(data, cols, "row")) .+ 1
    colidx = Int.(_grid_column(data, cols, "col")) .+ 1
    vals = _grid_column(data, cols, "value")
    (maximum(rows) <= n_grid && maximum(colidx) <= n_grid) ||
        error("cholesky indices exceed the $(n_grid)-point MESMER land grid")
    return sparse(rows, colidx, vals, n_grid, n_grid)
end

# 【grid_fingerprint】几何坐标的顺序敏感指纹（open_issues.md R20）。
#
# 必须与 03Claim/scripts/grid_fingerprint.py 逐位一致，因为指纹由 Python 侧的
# build_mesmer_country_weights.py 写下、由这里校验。可复现性靠两条约定保证：
#   ① 坐标先化成**整数微度**再哈希——float64 的十进制舍入在两种语言里未必同解，
#      整数运算则逐位一致；
#   ② 显式小端 int64、lat/lon 交错成 (n,2)，两侧按同样字节序拼装。
# 顺序敏感是刻意的：land_idx 就是行号，对调两行就是另一套网格。
const GRID_FINGERPRINT_SCHEMA = "mesmer-grid-v1"
const MICRODEGREE = 1_000_000

function grid_fingerprint(lat::AbstractVector, lon::AbstractVector)
    length(lat) == length(lon) ||
        error("lat/lon length mismatch: $(length(lat)) vs $(length(lon))")
    n = length(lat)
    buf = IOBuffer()
    write(buf, Vector{UInt8}("$(GRID_FINGERPRINT_SCHEMA) n=$(n)\n"))
    for i in 1:n
        la = Float64(lat[i])
        lo = Float64(lon[i])
        isfinite(la) && isfinite(lo) || error("coordinates contain NaN or Inf at row $i")
        abs(la) <= 90.0 + 1e-9 || error("latitude outside [-90, 90] at row $i")
        lo = mod(lo + 180.0, 360.0) - 180.0
        ila = round(Int64, la * MICRODEGREE)
        ilo = round(Int64, lo * MICRODEGREE)
        ilo == 180 * MICRODEGREE && (ilo = -180 * MICRODEGREE)
        write(buf, htol(ila))
        write(buf, htol(ilo))
    end
    return bytes2hex(sha256(take!(buf)))[1:16]
end

# 【geometry_fingerprint】对 climate_grid_geometry*.csv 按文件行序取指纹。
function geometry_fingerprint(path::AbstractString)
    isfile(path) || error("missing grid geometry file: $path")
    data, cols = _read_numeric_csv(path)
    haskey(cols, "lat") && haskey(cols, "lon") ||
        error("$path has no lat/lon columns")
    return grid_fingerprint(view(data, :, cols["lat"]), view(data, :, cols["lon"]))
end

# 【check_grid_weights_fingerprint】W 是否按当前几何算出来的（R20）。
#
# 原有守卫只查 land_idx 是否越界，那挡不住"同维度、坐标错位"：改一下岛屿补格
# 列表而格点总数凑巧不变时，每个 land_idx 指向的地点全变了，旧 W 照样通过，
# 于是每个国家的权重被投到错误的格点上——跑得出数、不报错、结果全错。
#
# 缺 manifest 与指纹不符同样按陈旧处理：两者的后果一样，都无法证明对齐。
function check_grid_weights_fingerprint(unshared_dir::AbstractString)
    manifest = joinpath(unshared_dir, "climate_grid_weights_manifest.csv")
    geometry = joinpath(unshared_dir, "climate_grid_geometry.csv")
    # 几何缺失必须报错，不能放行。这里原本写的是 `isfile(geometry) || return nothing`
    # 并注了一句"几何缺失由上游加载器各自报错"——**那句话是错的**：这条路径上没有
    # 任何别处读几何，于是删掉 climate_grid_geometry.csv 就能整个绕过指纹检查，
    # 维度正确的旧权重照常加载。这正是 R8/R17 反复记录的"守卫存在但不能失败"，
    # 由外部审查在 2026-08-06 查出（codex/gpt-5.6-sol）。
    isfile(geometry) || error(
        "climate_grid_geometry.csv not found in $unshared_dir. The weight matrix " *
        "cannot be shown to match the current grid without it, and land_idx is " *
        "positional — a mismatch is silent, not an error (open_issues.md R20). " *
        "Regenerate the grid parameters with data/mesmer_data/export_grid_params.py.")
    isfile(manifest) || error(
        "climate_grid_weights.csv has no fingerprint manifest at $manifest. " *
        "It predates the R20 guard and cannot be shown to match the current grid; " *
        "rebuild it with scripts/build_mesmer_country_weights.py (open_issues.md R20).")
    rows = _read_csv_dicts(manifest)
    recorded = ""
    for row in rows
        row["key"] == "grid_fingerprint" && (recorded = row["value"])
    end
    want = geometry_fingerprint(geometry)
    recorded == want || error(
        "stale climate_grid_weights.csv: its manifest records grid_fingerprint=" *
        "$(repr(recorded)) but climate_grid_geometry.csv now fingerprints as $(repr(want)). " *
        "Matching row counts are not enough — land_idx is positional, so shifted " *
        "coordinates would silently mis-project every country (open_issues.md R20). " *
        "Rebuild with scripts/build_mesmer_country_weights.py.")
    return want
end

# 【load_grid_weights】人口权重矩阵 W（regions × land grid）。逐国权重和必须为 1，
# 且每国至少一个正权重——旧链路正是在这里静默补零的。
function load_grid_weights(unshared_dir::AbstractString, regions, n_grid::Int)
    check_grid_weights_fingerprint(unshared_dir)
    rows = _read_csv_dicts(joinpath(unshared_dir, "climate_grid_weights.csv"))
    position = Dict(string(r) => i for (i, r) in enumerate(regions))
    ri, ci, vv = Int[], Int[], Float64[]
    for row in rows
        haskey(position, row["region"]) || continue
        push!(ri, position[row["region"]])
        push!(ci, parse(Int, row["land_idx"]) + 1)
        push!(vv, parse(Float64, row["weight"]))
    end
    W = sparse(ri, ci, vv, length(regions), n_grid)
    sums = vec(sum(W, dims=2))
    for (i, r) in enumerate(regions)
        isapprox(sums[i], 1.0; atol=1e-8) ||
            error("weight row for region $r sums to $(sums[i]); expected 1")
    end
    return W
end

# 【load_mesmer_grid】一次性加载全部网格参数与权重矩阵。
function load_mesmer_grid(data_dir::AbstractString, regions)
    unshared = joinpath(data_dir, "unshared_parameters")
    lt = load_grid_lt(unshared)
    lv = load_grid_lv(unshared)
    n_grid = length(lt.beta)
    length(lv.coef_gvtas) == n_grid ||
        error("climate_grid_lt.csv and climate_grid_lv.csv disagree on the grid size")
    return (n_grid=n_grid, lt=lt, lv=lv, gv=load_grid_gv(unshared),
            chol=load_innov_cholesky(unshared, n_grid),
            W=load_grid_weights(unshared, regions, n_grid))
end

# ----------------------------------------------------------------------------
# 二、聚合与模拟
# ----------------------------------------------------------------------------

# 【national_lt_params】把网格 β/γ 人口加权投影到国家层，得到 climateregional
# 组件需要的 region_lt_coef / region_lt_inter。
function national_lt_params(grid)
    return (ltcoef=Vector(grid.W * grid.lt.beta), ltinter=Vector(grid.W * grid.lt.gamma))
end

# 【_draw_ar_scalar】复刻 mesmer 的标量 AR(p) 抽样（全球变率）。上游递推见
# mesmer/stats/_auto_regression.py:494-501：前 ar_order+1 个时点保持为 0，
# 从 t = ar_order+1（0-based）起递推，最后丢弃 buffer 段。
function _draw_ar_scalar(order, intercept, coefs, std_innovs, n_ts, rng;
                         buffer::Int=MESMER_GV_BUFFER)
    total = n_ts + buffer
    out = zeros(Float64, total)
    innov = randn(rng, total) .* std_innovs
    for t in (order + 2):total          # 0-based range(order+1, total) -> 1-based
        acc = intercept
        for (lag, c) in enumerate(coefs)
            acc += c * out[t - lag]
        end
        out[t] = acc + innov[t]
    end
    return out[(buffer + 1):end]
end

# 【_draw_ar1_grid】复刻 mesmer 的网格 AR(1) + 空间相关新息抽样。
# 新息 ε_t = L z_t（z 标准正态），L 为 loc_ecov_AR1_innovs 的 Cholesky 因子。
# 与上游一致：前两个时点保持为 0，buffer 段丢弃。
function _draw_ar1_grid(ar1_int, ar1_coef, chol, n_ts, rng;
                        buffer::Int=MESMER_LV_BUFFER)
    n_grid = length(ar1_int)
    total = n_ts + buffer
    out = Matrix{Float64}(undef, n_ts, n_grid)
    prev = zeros(Float64, n_grid)
    z = Vector{Float64}(undef, n_grid)
    for t in 1:total
        if t > 2
            randn!(rng, z)
            prev = ar1_int .+ ar1_coef .* prev .+ chol * z
        end
        if t > buffer
            @inbounds out[t - buffer, :] = prev
        end
    end
    return out
end

# 【simulate_grid_lv】网格级局地变率完整分解：
#   ν_{i,t} = coef_gvtas_i · gv_t + a_{i,t}
# gv 为 AR(p) 全球变率，a 为 AR(1)+空间相关新息过程。参见 create_emus_lv.py 的
# OLS_AR1_sci 组合；params_lv 无 intercept 键（fit_intercept=False），故 OLS
# 分量不含截距。
function simulate_grid_lv(grid, n_ts::Int, rng)
    c = simulate_grid_lv_components(grid, n_ts, rng)
    nu = c.loc
    @inbounds for t in 1:n_ts
        nu[t, :] .+= grid.lv.coef_gvtas .* c.gv[t]
    end
    return nu
end

# 【simulate_grid_lv_components】把上式的两项**分开**返回：
#   gv  —— 全球年际变率，标量时间序列（ENSO 一类，各国同向）
#   loc —— 逐格点 AR(1) + 空间相关新息，纯特异性局地天气
# 完整 ν = loc .+ coef_gvtas' .* gv，与 simulate_grid_lv 逐值相同。
#
# 为何需要分开（open_issues.md P1-19）：回归端 02_lp_regression.do §6.5 把局地
# 温度对 g_shock（全球温度的 Hamilton 残差 = 全球变率）做了逐国斜率的正交化，
# 即"全球变率驱动的那部分局地温度"其损害效应被归入**全球**渠道 ζ_g。若把完整
# ν 喂给 ζ_l，这部分就用了错误渠道的系数——实测占国家级 ν 方差的中位 19%。
# 故损害组件消费 loc，而 gv 单独并入全球渠道的驱动项。
#
# **抽取顺序必须与合并版严格一致**（先 gv 后网格 AR(1)），否则同一 seed 会得到
# 不同实现，与 climate_region_lv.csv 的生产基准对不上。
function simulate_grid_lv_components(grid, n_ts::Int, rng)
    gv = _draw_ar_scalar(grid.gv.order, grid.gv.ar_int, grid.gv.coefs,
                         grid.gv.ar_std, n_ts, rng)
    loc = _draw_ar1_grid(grid.lv.ar1_int, grid.lv.ar1_coef, grid.chol, n_ts, rng)
    return (gv=gv, loc=loc)
end

# 【simulate_region_lv】生成一条国家级局地天气变率路径（time × regions）。
# 历史段（< variability_start_year）置零，与旧口径一致。
# seed 必须显式给出：SCC 的 Base/Pulse 两轨要求共用同一条实现，且结果需可复现。
function simulate_region_lv(grid, years; seed::Integer,
                            variability_start_year::Int=MESMER_VARIABILITY_START_YEAR)
    return simulate_region_lv_components(grid, years; seed=seed,
        variability_start_year=variability_start_year).total
end

# 【simulate_region_lv_components】一次抽样，返回国家级的三条东西：
#   total —— 完整 ν，喂给 climateregional 生成国别温度产出（口径不变）
#   loc   —— 纯特异性局地天气，喂给损害组件的局地渠道 ζ_l
#   gv    —— 全球年际变率标量，并入损害组件全球渠道 ζ_g 的驱动项
# 三者满足 total = loc + W·(coef_gvtas .* gv)，见 P1-19。
#
# gv 与 ν 同样在 variability_start_year 之前置零：历史段的国别温度直接取观测
# （P1-16），此时若全球渠道另外加一份变率，就与观测段自相矛盾。
function simulate_region_lv_components(grid, years; seed::Integer,
                                       variability_start_year::Int=MESMER_VARIABILITY_START_YEAR)
    years = collect(years)
    nR = size(grid.W, 1)
    future = findall(y -> y >= variability_start_year, years)
    total = zeros(Float64, length(years), nR)
    loc = zeros(Float64, length(years), nR)
    gv = zeros(Float64, length(years))
    isempty(future) && return (total=total, loc=loc, gv=gv)
    rng = MersenneTwister(seed)
    c = simulate_grid_lv_components(grid, length(future), rng)
    Wt = transpose(grid.W)
    @inbounds loc[future, :] = c.loc * Wt
    @inbounds gv[future] = c.gv
    gvgrid = c.gv * transpose(grid.lv.coef_gvtas)
    @inbounds total[future, :] = loc[future, :] .+ gvgrid * Wt
    return (total=total, loc=loc, gv=gv)
end

# 【write_region_lv】把一条国家级变率实现写成长表 CSV，供溯源与出图使用。
# 调用方应把 seed 写进文件名，杜绝旧链路"不知道用了哪个 emu"的溯源断裂。
function write_region_lv(path::AbstractString, years, regions, lv::AbstractMatrix)
    open(path, "w") do io
        println(io, "year,region,value")
        for (ti, y) in enumerate(years), (ri, r) in enumerate(regions)
            println(io, "$(y),$(r),$(lv[ti, ri])")
        end
    end
    return path
end

# 【_MESMER_GRID_CACHE】按 data_dir 缓存已加载的网格参数。Cholesky 因子约
# 110 MB / 400 万非零，嵌套 MC 每次 get_model 都重读会非常浪费。
const _MESMER_GRID_CACHE = Dict{String,Any}()

function cached_mesmer_grid(data_dir::AbstractString, regions)
    key = string(abspath(data_dir), "|", join(regions, ","))
    return get!(_MESMER_GRID_CACHE, key) do
        load_mesmer_grid(data_dir, regions)
    end
end

# 【regenerate_national_mesmer_inputs】由网格参数 + 人口权重矩阵重新生成
# climateregional 组件消费的三个国家级 CSV：
#   climate_ltcoef.csv / climate_ltinter.csv —— W·β 与 W·γ，确定性，无 seed；
#   climate_region_lv.csv                    —— 一条 seeded 局地变率实现。
# 这样模型端的加载路径完全不变，而口径从"面积加权 + RegionID 中间层"换成
# "人口加权 + 网格直投"。seed 会一并写进 climate_region_lv_manifest.csv，
# 修掉旧链路"不知道用了哪个 emu"的溯源断裂。
#
# 注意：SCC 的 Base/Pulse 两轨必须共用同一条变率实现（指南 §3.6 局地渠道
# 在边际中相消的前提）。因此本函数应在 get_model 之前调用一次，不可在
# Base 与 Pulse 之间重新生成。
#
# output_dir 默认写回 data_dir/unshared_parameters；做新旧差异审计时先指向
# 临时目录，确认无误后再正式覆盖（方案第 7 条：新文件验证通过前不得覆盖
# 现有国家 climate 文件）。
function regenerate_national_mesmer_inputs(data_dir::AbstractString, regions, years;
                                           seed::Integer,
                                           variability_start_year::Int=MESMER_VARIABILITY_START_YEAR,
                                           output_dir::AbstractString="")
    grid = cached_mesmer_grid(data_dir, regions)
    unshared = isempty(output_dir) ? joinpath(data_dir, "unshared_parameters") : String(output_dir)
    mkpath(unshared)
    lt = national_lt_params(grid)

    open(joinpath(unshared, "climate_ltcoef.csv"), "w") do io
        println(io, "region,value")
        for (i, r) in enumerate(regions)
            println(io, "$(r),$(lt.ltcoef[i])")
        end
    end
    open(joinpath(unshared, "climate_ltinter.csv"), "w") do io
        println(io, "region,value")
        for (i, r) in enumerate(regions)
            println(io, "$(r),$(lt.ltinter[i])")
        end
    end

    # 一次抽样、三个产物，必须配套使用（open_issues.md P1-19）：
    #   climate_region_lv.csv       完整 ν → climateregional，决定国别温度产出
    #   climate_region_lv_local.csv 纯局地天气 → 损害组件的局地渠道 ζ_l
    #   climate_global_lv.csv       全球年际变率 → 损害组件全球渠道 ζ_g 的驱动项
    # 三者来自同一 seed 的同一次抽取，混用不同批次会破坏 total = loc + W·(c^gv·gv)。
    c = simulate_region_lv_components(grid, years; seed=seed,
                                      variability_start_year=variability_start_year)
    lv = c.total
    write_region_lv(joinpath(unshared, "climate_region_lv.csv"), years, regions, lv)
    write_region_lv(joinpath(unshared, "climate_region_lv_local.csv"), years, regions, c.loc)
    open(joinpath(unshared, "climate_global_lv.csv"), "w") do io
        println(io, "year,value")
        for (ti, y) in enumerate(years)
            println(io, "$(y),$(c.gv[ti])")
        end
    end

    open(joinpath(unshared, "climate_region_lv_manifest.csv"), "w") do io
        println(io, "key,value")
        println(io, "seed,$(seed)")
        println(io, "first_year,$(first(years))")
        println(io, "last_year,$(last(years))")
        println(io, "variability_start_year,$(variability_start_year)")
        println(io, "n_regions,$(length(regions))")
        println(io, "n_land_gridpoints,$(grid.n_grid)")
        println(io, "gv_ar_order,$(grid.gv.order)")
        println(io, "rng,MersenneTwister")
    end
    return (ltcoef=lt.ltcoef, ltinter=lt.ltinter, region_lv=lv,
            region_lv_local=c.loc, global_lv=c.gv)
end

# 【grid_warming】给定全球温度路径，还原网格级升温的强迫响应项
# β_i (T_global − T_ref) + γ_i，供 MimiCLAIM.ipynb 出网格化温升图。
# 不含天气噪声：论文图通常画强迫响应。
function grid_warming(grid, T_global, T_global_ref)
    return [grid.lt.beta[i] * (T_global[t] - T_global_ref) + grid.lt.gamma[i]
            for t in eachindex(T_global), i in 1:grid.n_grid]
end

# ----------------------------------------------------------------------------
# 三、存档：旧 RegionID 聚合链路的读取器
# ----------------------------------------------------------------------------
# 以下两个函数读取 data/mesmer_data/mesmer_results/ 下按 2,388 个次国家
# RegionID 聚合的历史文件。它们是当前 SCC 基准的唯一来源，保留下来只为做
# 新旧对比与解释结果差异，**不参与模型运行**（见 open_issues.md P1-3）。

function load_region_lt(esm)
    path = "data/mesmer_data/mesmer_results/"
    region_lt_coef  = DataFrame(load(joinpath(@__DIR__, "../"*path*esm*"_regional_T_coef.csv")))
    region_lt_inter = DataFrame(load(joinpath(@__DIR__, "../"*path*esm*"_regional_T_inter.csv")))
    return region_lt_coef, region_lt_inter
end

function load_region_lv(esm, emu)
    path = "data/mesmer_data/mesmer_results/"
    region_lv = DataFrame(load(joinpath(@__DIR__, "../"*path*esm*"_regional_variability_"*string(emu)*".csv")))
    return region_lv
end
