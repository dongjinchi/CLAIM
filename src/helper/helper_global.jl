# ============================================================================
# 【helper_global.jl】全局层参数与装配 helper——负责损害模块以外的所有
# "模型级"事务：通用 CSV 读取器、维度/区域设置、组件装配与跨组件连线、
# shared_parameters 注册，以及 socioeconomic / climateregional 两个基线
# 组件的 unshared 参数加载。
# 损害模块（components_damage/claim_damage）的参数由 helper_damage.jl 负责；
# FaIR 气候链参数由 helper_fairv2.jl 负责。
# 情景数据说明：get_model 的 scenario（ssp119/ssp126/ssp245/ssp370/ssp585）
# 经 _socio_ssp_from_scenario 映射到 IIASA SSP v3.2 社会经济文件
# population_sspX.csv / socioeconomic_ypcgrowth_sspX.csv（PWT 11.0 历史
# 1950–2024 + 情景 2025–2500，由 01data/ssp/ 生成）。_load_time_region_matrix
# 对超出 CSV 最后年份的模型年份自动"平推"最后一个观测值。
# ============================================================================

const CLIMATE_BASE_YEAR = 1950
const T_GLOBAL_REF_1950 = 0.2749925941187845

# 【CLIMATE_OBSERVED_END_YEAR】历史观测段的末年。1950–2024 的国家温度直接取
# Berkeley 人口加权观测（与 LP 的 lctmp_bkly_pw 同源），2025 起才由 MESMER
# 模拟。选 2024 是因为它同时是 variability_start_year 与 bau_ref_year 的前一年，
# 且 Berkeley 的 2024 是完整年份（见 open_issues.md P1-16）。
const CLIMATE_OBSERVED_END_YEAR = 2024

# 【CLIMATE_ANCHOR_WINDOW】前向段水平锚的观测窗口长度（年）。取 30 与
# T_state 的气候态窗口、LP 的 barT_s 窗口一致；用单年值会把拼接年的天气异常
# 永久固化进前向趋势线。
const CLIMATE_ANCHOR_WINDOW = 30

function climate_global_ref(; climate_base_year::Int=CLIMATE_BASE_YEAR)
    climate_base_year == CLIMATE_BASE_YEAR || error("Only climate_base_year=1950 is currently supported")
    return T_GLOBAL_REF_1950
end

# ----------------------------------------------------------------------------
# 一、通用 CSV 读取器
# ----------------------------------------------------------------------------

# 【_read_csv_dicts】读取简单 CSV 文件为 Dict 行向量；B&K 参数表均采用
# 无引号逗号分隔的数值/字符串列，避免引入 JSON 或额外依赖。
function _read_csv_dicts(path::AbstractString)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    rows = Vector{Dict{String,String}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = split(strip(line), ',')
        push!(rows, Dict(header[i] => values[i] for i in eachindex(header)))
    end
    return rows
end

function _load_region_vector(path::AbstractString, regions)
    rows = _read_csv_dicts(path)
    values = Dict(row["region"] => parse(Float64, get(row, "value", get(row, "ypc0", get(row, "tempbase", get(row, "temp1990", "0"))))) for row in rows)
    return [values[string(r)] for r in regions]
end

# 【_load_time_region_matrix】读取 time × region 输入（open_issues.md R19-b）。
#
# 这个函数原本对**任何**残缺输入都静默补值，第三轮外部审查（2026-08-06）点名它
# 是 R19 同类缺陷的主体：
#   · 整国缺失 → 整列写 `default`（0.0）。从 `population_ssp2.csv` 删掉 USA，
#     模型把 USA 全期人口设为零，跑完不报错；
#   · 中间缺年 / 尾部缺年 → LOCF（最后观测值前推）。删掉某国 2030 年的增长率，
#     它会静默沿用 2029 年的值。
#
# 但**前导缺口是正当的**，不能一并禁掉：生产 `population_ssp2.csv` 里 130 个国家
# 缺 1950–1969（ARM/AZE 缺到 1989——苏联解体前不作为独立经济体存在），PWT 本就
# 没有那些年份。对前导缺口用首值回填是有意的。
#
# 故判据按缺口位置区分，而不是"有缺就报错"：
#   ① 整国缺失            → 一律报错（没有哪种情形下"这国不存在"是对的）
#   ② 同一国的年份不连续  → 报错（中间缺口）
#   ③ 序列末年 < 模型末年 → 报错（尾部缺口，即 LOCF 外推）
#   ④ 序列首年 > 模型首年 → 放行，用首值回填（前导缺口）
#   ⑤ (国, 年) 重复       → 报错
#
# `strict=false` 只留给合成数据测试，且必须显式点名。
# `required_through` 让调用方声明"这份文件只需覆盖到哪一年"。默认是模型末年；
# 观测温度这类**刻意只覆盖历史段**的输入把它设成拼接年，于是尾部判据仍然生效，
# 只是标尺换成了正确的那一把——而不是靠 `strict=false` 把三条判据一起关掉。
function _load_time_region_matrix(path::AbstractString, years, regions;
                                  default::Float64=0.0, strict::Bool=true,
                                  required_through::Union{Int,Nothing}=nothing)
    rows = _read_csv_dicts(path)
    raw = Dict{String,Vector{Tuple{Int,Float64}}}()
    for row in rows
        region = row["region"]
        push!(get!(raw, region, Tuple{Int,Float64}[]), (parse(Int, row["year"]), parse(Float64, row["value"])))
    end
    fname = basename(path)
    last_year = required_through === nothing ? maximum(years) : required_through
    out = Array{Float64}(undef, length(years), length(regions))
    for (ri, r) in enumerate(regions)
        series = sort(get(raw, string(r), Tuple{Int,Float64}[]), by=x -> x[1])
        if isempty(series)
            strict && error(
                "$fname has no rows for region $r. A missing country would be filled " *
                "with $default for every year and the run would still succeed " *
                "(open_issues.md R19-b). Regenerate the file, or pass strict=false " *
                "if this is a synthetic-data test.")
            out[:, ri] .= default
            continue
        end
        if strict
            ys = first.(series)
            length(unique(ys)) == length(ys) || error(
                "$fname has duplicate (region, year) rows for $r; the later row would " *
                "silently win (open_issues.md R19-b).")
            ys == collect(ys[1]:ys[end]) || error(
                "$fname is missing interior years for region $r (has $(ys[1])–$(ys[end]) " *
                "with $(ys[end] - ys[1] + 1 - length(ys)) gaps). Gaps are carried forward " *
                "from the last observation, which looks like data (open_issues.md R19-b).")
            ys[end] >= last_year || error(
                "$fname ends at $(ys[end]) for region $r but the model runs to $last_year. " *
                "The tail would be filled by carrying $(ys[end]) forward for " *
                "$(last_year - ys[end]) years (open_issues.md R19-b).")
        end
        for (ti, y) in enumerate(years)
            idx = findlast(x -> x[1] <= y, series)
            out[ti, ri] = idx === nothing ? series[1][2] : series[idx][2]
        end
    end
    return out
end

function _load_time_vector(path::AbstractString, years)
    rows = _read_csv_dicts(path)
    series = sort([(parse(Int, row["year"]), parse(Float64, row["value"])) for row in rows], by=x -> x[1])
    return [begin
        idx = findlast(x -> x[1] <= y, series)
        idx === nothing ? series[1][2] : series[idx][2]
    end for y in years]
end

function _parse_float(row, key)
    return parse(Float64, row[key])
end

function _parse_int(row, key)
    return parse(Int, row[key])
end

# ----------------------------------------------------------------------------
# 二、维度、区域与组件装配
# ----------------------------------------------------------------------------

function load_national_regions(shared_dir::AbstractString=joinpath(@__DIR__, "..", "..", "data", "shared_parameters"))
    rows = _read_csv_dicts(joinpath(shared_dir, "regions_national.csv"))
    return [row["region"] for row in rows]
end

function setup_claim_dimensions!(m, years, regions)
    set_dimension!(m, :time, years)
    set_dimension!(m, :regions, regions)
    set_dimension!(m, :horizon_dim, 0:10)
end

function active_damage_component(damage_function::Symbol)
    damage_function == :BHM_global && return :bhm_damage
    damage_function in (:BK_global, :BK_region) && return :claim_damage
    error("damage_function must be :BHM_global, :BK_global, or :BK_region")
end

function validate_damage_controls(damage_function::Symbol, income_state::Symbol, capital_response::Bool)
    damage_function in (:BHM_global, :BK_global, :BK_region) || error("damage_function must be :BHM_global, :BK_global, or :BK_region")
    income_state in (:none, :dynamic) || error("income_state must be :none or :dynamic")
    if damage_function in (:BHM_global, :BK_global) && income_state != :none
        error("damage_function=$damage_function requires income_state=:none")
    end
    if damage_function == :BHM_global && capital_response
        error("damage_function=:BHM_global currently requires capital_response=false")
    end
    return nothing
end

function add_claim_components!(m; damage_component::Symbol=:claim_damage)
    add_comp!(m, socioeconomic)
    add_comp!(m, climateregional)
    if damage_component == :claim_damage
        add_comp!(m, claim_damage)
    elseif damage_component == :bhm_damage
        add_comp!(m, bhm_damage)
    else
        error("unsupported damage_component=$damage_component")
    end
end

# 【connect_claim_components!】跨组件连线（climateregional→损害组件、
# socioeconomic→损害组件、shared population→各消费组件）。
function connect_claim_components!(m; damage_component::Symbol=:claim_damage)
    connect_param!(m, damage_component, :regtmp, :climateregional, :regtmp)
    connect_param!(m, damage_component, :ypc_baseline, :socioeconomic, :ypc)
    connect_param!(m, :socioeconomic, :population, :population)
    connect_param!(m, damage_component, :population, :population)
end

# ----------------------------------------------------------------------------
# 三、shared_parameters 注册
# ----------------------------------------------------------------------------

# 【_socio_ssp_from_scenario】排放情景 → 社会经济情景文件名后缀。
# 与 data/fair_data 中五种 RCMIP 情景一一对应：首位数字即 SSP 编号
#（ssp119/ssp126→ssp1，ssp245→ssp2，ssp370→ssp3，ssp585→ssp5）。
function _socio_ssp_from_scenario(scenario::AbstractString)
    return "ssp" * scenario[4]
end

# 【_socio_scenario_file】按社会经济情景后缀解析文件名。
function _socio_scenario_file(base::AbstractString, ssp::AbstractString)
    return base * "_" * lowercase(ssp) * ".csv"
end

function add_claim_shared_params!(m; data_dir::AbstractString, years, regions, ssp::AbstractString)
    shared_dir = joinpath(data_dir, "shared_parameters")
    popfile = _socio_scenario_file("population", ssp)
    add_shared_param!(m, :population, _load_time_region_matrix(joinpath(shared_dir, popfile), years, regions), dims=[:time, :regions])
end

# ----------------------------------------------------------------------------
# 四、全局层 unshared 参数加载（socio / climateregional 分节）
# ----------------------------------------------------------------------------

# 【update_global_params!】加载 socioeconomic 与 climateregional 两个基线
# 组件的参数。climate=:csv 时 inputtemp 由 shared_parameters/climate_T_global.csv
# 外生给定；climate=:fair 时 inputtemp 已在 get_model 中连接到 temperature:T，
# 此处不覆盖。
function update_global_params!(m; data_dir::AbstractString, years, regions, climate::Symbol=:csv,
                               ssp::AbstractString, start_year::Int=CLIMATE_BASE_YEAR,
                               allow_legacy_climate_anchor::Bool=false)
    shared_dir = joinpath(data_dir, "shared_parameters")
    unshared_dir = joinpath(data_dir, "unshared_parameters")

    # ---- socioeconomic：初始人均收入 + 外生增长率路径（可选 SSP 情景）----
    update_param!(m, :socioeconomic, :ypc0, _load_region_vector(joinpath(unshared_dir, "socioeconomic_ypc0.csv"), regions))
    growthfile = _socio_scenario_file("socioeconomic_ypcgrowth", ssp)
    update_param!(m, :socioeconomic, :ypcgrowth, _load_time_region_matrix(joinpath(unshared_dir, growthfile), years, regions))

    # ---- climateregional：MESMER 降尺度系数与区域天气变率 ----
    update_param!(m, :climateregional, :region_lt_coef, _load_region_vector(joinpath(unshared_dir, "climate_ltcoef.csv"), regions))
    update_param!(m, :climateregional, :region_lt_inter, _load_region_vector(joinpath(unshared_dir, "climate_ltinter.csv"), regions))
    update_param!(m, :climateregional, :region_lv, _load_time_region_matrix(joinpath(unshared_dir, "climate_region_lv.csv"), years, regions))
    update_param!(m, :climateregional, :tempbase, _load_region_vector(joinpath(unshared_dir, "climate_tempbase.csv"), regions))
    update_param!(m, :climateregional, :T_global_ref, climate_global_ref())

    # ---- 历史段的观测国家温度（open_issues.md P1-16 / R19）----
    # 1950–2024 直接用 Berkeley 人口加权观测，不再由 MESMER 重构；前向段的水平
    # 锚定在观测窗口均值上。
    #
    # 这里曾是 `isfile ? 观测 : @warn + legacy 重构`。R19 指出那是个**静默降级**：
    # P1-15/P1-16 之所以要修，正是因为 legacy 分支把水平锚死在 1950 单年上，而
    # `@warn` 在批量运行里滚过去没人看，跑出来的 SCC 与生产口径不是一回事却毫无
    # 标记。缺观测不是一种"配置"，是生产输入不完整，必须硬失败。
    #
    # 两种失效方式都要挡，因为它们导致同一个后果（regtmp 全程由模型重构）：
    #   ① 文件不存在；
    #   ② 文件存在但整个观测段落在模型时间轴之外（start_year > 拼接年）。
    #      ② 尤其阴险：`_load_time_region_matrix` 是 LOCF，会把 2024 年的值一路
    #      前推填满 2025+ 的每一行，看上去"有观测"，而组件那侧 `has_observed`
    #      为假，实际走的仍是 legacy。文件在、数字在、结果错。
    #
    # legacy 分支只保留给合成数据测试，且必须由调用方显式点名
    # `allow_legacy_climate_anchor=true`——名字里写明它是 legacy 锚，
    # 不是"没找到文件时的自动兜底"。
    observed_path = joinpath(unshared_dir, "climate_regtmp_observed.csv")
    has_file = isfile(observed_path)
    covers_horizon = start_year <= CLIMATE_OBSERVED_END_YEAR
    if has_file && covers_horizon
        # R19-c：锚窗被截断的配置必须挡住。组件用 `max(1, lo)` 截断锚窗，于是
        # start_year=2020 时"30 年锚"实际只有 2020–2024 五年，start_year=2024
        # 更退化成单年锚——**那正是 P1-15 判定为错误、P1-16 修掉的口径**，只是
        # 换了个入口悄悄回来。
        #
        # 判据只在锚**真的会被用到**时生效：模型末年不跨过拼接年时，前向段根本
        # 不存在，短窗无害（多数合成测试是 end_year=2021 这种）。
        anchor_start = CLIMATE_OBSERVED_END_YEAR - CLIMATE_ANCHOR_WINDOW + 1
        if maximum(years) > CLIMATE_OBSERVED_END_YEAR && start_year > anchor_start
            error("start_year=$start_year truncates the $(CLIMATE_ANCHOR_WINDOW)-year " *
                  "anchor window to $(CLIMATE_OBSERVED_END_YEAR - start_year + 1) years " *
                  "($(start_year)–$(CLIMATE_OBSERVED_END_YEAR)), and this run crosses the " *
                  "splice year so the anchor is actually used. A short window re-creates " *
                  "the single-year anchor that P1-15/P1-16 removed — 2024 alone sits " *
                  "+0.368 °C above trend (open_issues.md R19-c). Use start_year <= " *
                  "$anchor_start, or end_year <= $CLIMATE_OBSERVED_END_YEAR.")
        end
        # 这份文件**刻意**只覆盖 1950–2024：组件在拼接年之后不再读它。故尾部判据
        # 的标尺是"历史段与模型时间轴的交集"，既不是模型末年（那会要求它覆盖到
        # 2150），也不是拼接年（那会要求 end_year=2022 的短跑也备到 2024）。
        # 其余三条判据照常生效。
        observed = _load_time_region_matrix(
            observed_path, years, regions;
            required_through=min(CLIMATE_OBSERVED_END_YEAR, maximum(years)))
        update_param!(m, :climateregional, :regtmp_observed, observed)
        update_param!(m, :climateregional, :observed_end_year, CLIMATE_OBSERVED_END_YEAR)
    elseif allow_legacy_climate_anchor
        reason = has_file ?
            "start_year=$start_year is past the observed splice year $CLIMATE_OBSERVED_END_YEAR" :
            "climate_regtmp_observed.csv not found in $unshared_dir"
        @warn "using the legacy MESMER reconstruction for the historical segment " *
              "($reason); this is not the production caliber (open_issues.md P1-15/P1-16/R19)"
        update_param!(m, :climateregional, :regtmp_observed, zeros(length(years), length(regions)))
        update_param!(m, :climateregional, :observed_end_year, start_year - 1)
    elseif !has_file
        error("climate_regtmp_observed.csv not found in $unshared_dir. The historical " *
              "segment must come from observations, not from a MESMER reconstruction " *
              "(open_issues.md P1-16). Regenerate it by passing --temperature_cache to " *
              "scripts/prepare_national_inputs.py. Synthetic-data tests that genuinely " *
              "want the retired reconstruction must pass allow_legacy_climate_anchor=true.")
    else
        error("start_year=$start_year is later than the observed splice year " *
              "$CLIMATE_OBSERVED_END_YEAR, so climate_regtmp_observed.csv covers none of " *
              "the model horizon and every year would be reconstructed instead of observed " *
              "(open_issues.md R19). Use start_year <= $CLIMATE_OBSERVED_END_YEAR, or pass " *
              "allow_legacy_climate_anchor=true if the legacy anchor is genuinely intended.")
    end
    update_param!(m, :climateregional, :model_start_year, start_year)
    update_param!(m, :climateregional, :anchor_window, CLIMATE_ANCHOR_WINDOW)

    if climate == :csv
        update_param!(m, :climateregional, :inputtemp, _load_time_vector(joinpath(shared_dir, "climate_T_global.csv"), years))
    end
end
