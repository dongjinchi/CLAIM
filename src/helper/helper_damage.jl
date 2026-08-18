# ============================================================================
# 【helper_damage.jl】损害模块（components_damage/claim_damage）专用 helper——
# 只负责 B&K 损害组件的参数需求：damage_*.csv 缓存加载、闭式加权 OLS、
# θ 抽样展开，以及 update_claim_damage_params! 参数装配。
# 通用 CSV 读取器与全局层参数见 helper_global.jl（须先 include）。
# ============================================================================

# 【BKW_STATE_WINDOW】T_state 的回看窗口长度（年）。**不是自由参数**：它由
# LP 端 barT_s 的 1950–1979 三十年气候态窗口钉死，改这个值就等于让前向模型的
# 气候态语义偏离回归端识别的对象（open_issues.md P1-4）。
const BKW_STATE_WINDOW = 30

# 【load_state_centers】读取 state_centers.csv（单行 CSV），返回中心化常数
# barT_center（=Stata 第8部分的 barT_W）与 lny0_center（=lny0_W）。
# 前向模型的 T_state/Y_state 必须与 LP 回归端使用同一套中心化口径。
function load_state_centers(path::AbstractString)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    values = split(strip(lines[2]), ',')
    data = Dict(header[i] => parse(Float64, values[i]) for i in eachindex(header))
    return data
end

function load_damage_state_centers(unshared_dir::AbstractString=joinpath(@__DIR__, "..", "..", "data", "unshared_parameters"))
    return load_state_centers(joinpath(unshared_dir, "damage_state_centers.csv"))
end

function _sorted_horizons(rows)
    return sort(unique(_parse_int(row, "horizon") for row in rows))
end

function _values_by_horizon(rows, value_key)
    horizons = _sorted_horizons(rows)
    values = Dict(_parse_int(row, "horizon") => _parse_float(row, value_key) for row in rows)
    return [values[h] for h in horizons]
end

# 【closed_form_ols】闭式加权 OLS 反演：给定动态 θ 序列(h=0..10)与预计算的
# Solow 产出基底 Y1/Y2，解出损害核的两个线性系数 (A, C)。
#
# ŷ_h = A·Y1[h] + C·Y2[h] 对 (A,C) 线性，故这是一个 2×2 的加权正规方程：
#
#     [ Σ w Y1²   Σ w Y1Y2 ] [A]   [ Σ w Y1 θ ]
#     [ Σ w Y1Y2  Σ w Y2²  ] [C] = [ Σ w Y2 θ ]
#
# **A 与 C 联合求解，残差覆盖 h=0..10**（open_issues.md P1-17）。
# 此前 A 由 h=0 解析锁定（A = θ₀/Y1[0]，因 k̂₀=0 使 ŷ₀ 只含即期项）。该等式
# 本身正确，但 θ(0) 恰是整条 IRF 上最不可靠的点——全球渠道 0/185 国 |t|>1.96，
# SNR_level 中位 0.96。锁死后 C 被迫去补偿这个噪声点。实测四变体统一对原始 θ
# 评分：A 锁定 R²_w 中位 0.6829，A 自由 0.7594，**185/185 国改善**。
# （另测了 Barnichon & Brownlees 2019 的 Smooth LP：跨 horizon 平滑 0/185 国
#   改善，收益全部来自 A 的锚点，故不引入。）
#
# 与上游 03_invert_damage.py 的 4 参数 NLS 口径一致：那里冻结形状前也是让
# A 参与拟合，这里冻结形状后 A 仍与 C 一起解。
#
# C 施加**峰高**边界，不是箱界（open_issues.md P1-12，2026-08-04）：
# 核为 ζ_s = A·e^{-Bs} + C·s^ν·e^{-Ds}，峰高 |C|·(ν/D)^ν·e^{-ν} ≤ 0.37
# ⟺ |C| ≤ 0.37·e^ν·(D/ν)^ν。该上限随该国的 D 变化，故由调用方传入
# （`c_max`，来自 damage_structural_params.csv 的 C_max_g / C_max_l，
# 由 precompute_basis.py 预计算）。旧版写死 M1 ±0.05 / M2 ±0.15：那只在
# ν=1、D=0.05 处等价于 0.37 的峰高意图，D 一变就偏，ν 一变更偏。
#
# A 不设界（极寒国 θ₀>0 的升温红利照收；|A| 实测中位约 0.008）。
# 对应 coupling_plan.md §1.4；每次调用约 11 次点积 + 一个 2×2 解，仍是微秒级，
# 支撑"每国每年每 MC 抽样"的动态重反演。**换 ν 不改变这一点**：ν 只进 Y2 的
# 预计算，正规方程的结构不变，MC 每次抽样仍不需要重做非线性拟合。
function closed_form_ols(theta, Y1, Y2, w, c_max)
    s11 = 0.0; s12 = 0.0; s22 = 0.0; b1 = 0.0; b2 = 0.0
    n = min(length(theta), length(Y1), length(Y2), length(w), 11)
    for idx in 1:n
        wh = isfinite(w[idx]) ? w[idx] : 0.0
        wh <= 0.0 && continue
        y1 = Y1[idx]; y2 = Y2[idx]; th = theta[idx]
        s11 += wh * y1 * y1
        s12 += wh * y1 * y2
        s22 += wh * y2 * y2
        b1 += wh * y1 * th
        b2 += wh * y2 * th
    end
    det = s11 * s22 - s12 * s12
    C_hi = abs(float(c_max))
    isfinite(C_hi) || throw(ArgumentError(
        "closed_form_ols: c_max must be finite (got $(c_max)). It comes from " *
        "damage_structural_params.csv::C_max_g/C_max_l; regenerate the basis " *
        "cache with precompute_basis.py if that column is missing."))
    if abs(det) > 1e-18
        A = (s22 * b1 - s12 * b2) / det
        C = (s11 * b2 - s12 * b1) / det
        if abs(C) > C_hi
            # 【C 撞界：必须固定 C 重解 A】（open_issues.md §6.2 R18）
            #
            # 把无约束解的违约坐标"剪掉"**不是**有约束问题的解。正确解法：
            # 先把 A 剖出——对任意给定 C，最优 A(C) = (b1 − C·s12)/s11。代回目标
            # 函数得到关于 C 的一维凸二次函数，其无约束极小点恰是上面的 C。在
            # 区间 [−C_hi, C_hi] 上最小化一维凸二次函数，解就是 clamp(C)。
            #
            # 所以**被剪的 C 是对的，错的只有 A**：A 必须跟着 C 一起更新。
            # 实测差距可达一倍：theta=[0,1]、Y1=[1,1]、Y2=[0,1]、等权、C 界 0.05 时
            # 旧实现给 (0, 0.05) 加权 SSR 0.9025，本式给 (0.475, 0.05) SSR 0.45125。
            #
            # 当前生产参数下 0/370 条撞界，故本分支不改变已发布结果；但 C_hi =
            # 0.37·e^ν·(D/ν)^ν 是 D 的函数，一旦 D 进入抽样（R23）撞界就会真的
            # 发生，届时若不修就会系统性劣于真解且不报错。
            C = clamp(C, -C_hi, C_hi)
            s11 > 1e-18 && (A = (b1 - C * s12) / s11)
        end
        return A, C
    end
    # 退化情形（Y1 与 Y2 共线，或权重几乎全零）：退回只解 A，C 置零。
    A = s11 > 1e-18 ? b1 / s11 : (length(theta) > 0 ? theta[1] : 0.0)
    return A, 0.0
end

# 【load_basis_cache】加载 data/unshared_parameters 中的 B&K damage_*.csv，
# 返回 basis["ISO3"]["G"/"L"] → (model, Z1, Z2, Y1, Y2, alpha, delta, g)
# 的嵌套字典。由 precompute_basis.py 预先生成。
function load_basis_cache(unshared_dir::AbstractString)
    basis_rows = _read_csv_dicts(joinpath(unshared_dir, "damage_basis.csv"))
    struct_rows = _read_csv_dicts(joinpath(unshared_dir, "damage_structural_params.csv"))
    weight_rows = _read_csv_dicts(joinpath(unshared_dir, "damage_weights.csv"))
    cache = Dict{String,Any}()
    for row in struct_rows
        region = row["region"]
        cache[region] = Dict{String,Any}(
            "G" => Dict{String,Any}("model" => _parse_int(row, "model_g"), "C_max" => _parse_float(row, "C_max_g"), "alpha" => _parse_float(row, "alpha"), "delta" => _parse_float(row, "delta"), "g" => _parse_float(row, "g"), "rho" => _parse_float(row, "rho")),
            "L" => Dict{String,Any}("model" => _parse_int(row, "model_l"), "C_max" => _parse_float(row, "C_max_l"), "alpha" => _parse_float(row, "alpha"), "delta" => _parse_float(row, "delta"), "g" => _parse_float(row, "g"), "rho" => _parse_float(row, "rho")),
        )
    end
    for region in keys(cache)
        for channel in ("G", "L")
            channel_rows = [row for row in basis_rows if row["region"] == region && row["channel"] == channel]
            for field in ("Z1", "Z2", "Y1", "Y2", "Y1_nocap", "Y2_nocap")
                rows = [row for row in channel_rows if row["field"] == field]
                cache[region][channel][field] = _values_by_horizon(rows, "value")
            end
        end
        for (field, channel) in (("w_g", "G"), ("w_l", "L"))
            rows = [row for row in weight_rows if row["region"] == region && row["channel"] == channel]
            cache[region][field] = _values_by_horizon(rows, "weight")
        end
    end
    return cache
end

# 【load_theta_draws】加载 data/unshared_parameters 中的 damage_theta_*.csv，
# 包含 theta_base（点估计，11×6）与 theta_draws（n_draws×11×6）。
function load_theta_draws(unshared_dir::AbstractString)
    cols = ["theta_g0", "theta_g1", "theta_g2", "theta_l0", "theta_l1", "theta_l2"]
    base_rows = sort(_read_csv_dicts(joinpath(unshared_dir, "damage_theta_base.csv")), by=row -> _parse_int(row, "horizon"))
    base = [[_parse_float(row, col) for col in cols] for row in base_rows]
    draws_path = joinpath(unshared_dir, "damage_theta_draws.csv")
    draws = Any[]
    if isfile(draws_path)
        draw_rows = _read_csv_dicts(draws_path)
        draw_ids = sort(unique(_parse_int(row, "draw_id") for row in draw_rows))
        for draw_id in draw_ids
            rows = sort([row for row in draw_rows if _parse_int(row, "draw_id") == draw_id], by=row -> _parse_int(row, "horizon"))
            push!(draws, [[_parse_float(row, col) for col in cols] for row in rows])
        end
    end
    return Dict("theta_base" => base, "theta_draws" => draws)
end

# 【matrix_from_channel】从 basis 缓存中按国家顺序提取某渠道(G/L)某字段
# （如 Y1/Z1）的 horizon 数组，拼成 (n_regions × 11) 矩阵，
# 供 update_param! 直接喂给 claim_damage 的基底参数。
function matrix_from_channel(cache, regions, channel, field)
    out = Array{Float64}(undef, length(regions), length(cache[string(regions[1])][channel][field]))
    for (i, r) in enumerate(regions)
        out[i, :] .= Float64.(cache[string(r)][channel][field])
    end
    return out
end

# 【vector_from_basis】从 basis 缓存中按国家顺序提取某渠道某标量字段
# （如 alpha/delta/g/model），返回向量。
function vector_from_basis(cache, regions, channel, field)
    return [cache[string(r)][channel][field] for r in regions]
end

# 【weights_from_cache】从 basis 缓存中按国家顺序提取逆方差权重数组
# （field 为 "w_g" 或 "w_l"，各 11 个 horizon），拼成 (n_regions × 11) 矩阵，
# 供 closed_form_ols 的加权投影使用。
function weights_from_cache(cache, regions, field)
    out = Array{Float64}(undef, length(regions), length(cache[string(regions[1])][field]))
    for (i, r) in enumerate(regions)
        out[i, :] .= Float64.(cache[string(r)][field])
    end
    return out
end

function claim_damage_controls(damage_function::Symbol, income_state::Symbol, capital_response::Bool)
    damage_function in (:BK_global, :BK_region) || error("B&K damage controls require damage_function=:BK_global or :BK_region")
    if damage_function == :BK_global && income_state != :none
        error("damage_function=:BK_global requires income_state=:none")
    end
    income_state in (:none, :dynamic) || error("income_state must be :none or :dynamic")
    return Dict(
        :temperature_effect => damage_function == :BK_region ? 1.0 : 0.0,
        :income_effect => income_state == :dynamic ? 1.0 : 0.0,
        :capital_effect => capital_response ? 1.0 : 0.0,
    )
end

function load_bhm_quadratic_params(unshared_dir::AbstractString)
    rows = _read_csv_dicts(joinpath(unshared_dir, "bhm_quadratic_params.csv"))
    values = Dict(row["parameter"] => parse(Float64, row["value"]) for row in rows)
    return values
end

function update_bhm_damage_params!(m; data_dir::AbstractString, years, regions, start_year::Int, bau_ref_year::Int)
    unshared_dir = joinpath(data_dir, "unshared_parameters")
    params = load_bhm_quadratic_params(unshared_dir)
    update_param!(m, :bhm_damage, :beta1, params["beta1"])
    update_param!(m, :bhm_damage, :beta2, params["beta2"])
    update_param!(m, :bhm_damage, :model_start_year, start_year)
    update_param!(m, :bhm_damage, :damage_start_year, bau_ref_year)
    update_param!(m, :bhm_damage, :region_lv, _load_time_region_matrix(joinpath(unshared_dir, "climate_region_lv.csv"), years, regions))
    update_param!(m, :bhm_damage, :tempbase, _load_region_vector(joinpath(unshared_dir, "climate_tempbase.csv"), regions))
end

# 【theta_parameters】把第 theta_id 次 pooled θ 抽样（theta_id=0 表示点估计）
# 展开为 claim_damage 需要的 6 个参数矩阵（theta_g0/g1/g2/l0/l1/l2，
# 各为 n_regions × 11）。注意 MC 扰动的是 pooled 系数——所有国家共享同一
# 抽样值，国家异质性由运行时 T_state/Y_state 的线性投影内生产生。
function theta_parameters(theta_payload, theta_id::Int, regions)
    base = theta_payload["theta_base"]
    theta = theta_id == 0 ? base : theta_payload["theta_draws"][theta_id]
    arr = Array{Float64}(undef, length(theta), 6)
    for h in eachindex(theta)
        arr[h, :] .= Float64.(theta[h])
    end
    return Dict(
        :theta_g0 => repeat(arr[:, 1]', length(regions), 1),
        :theta_g1 => repeat(arr[:, 2]', length(regions), 1),
        :theta_g2 => repeat(arr[:, 3]', length(regions), 1),
        :theta_l0 => repeat(arr[:, 4]', length(regions), 1),
        :theta_l1 => repeat(arr[:, 5]', length(regions), 1),
        :theta_l2 => repeat(arr[:, 6]', length(regions), 1),
    )
end

# 【update_claim_damage_params!】claim_damage 组件参数总装配。
# climate=:csv 时，T_global 由 shared_parameters/climate_T_global.csv 外生
# 给定，T0_global 取 bau_ref_year 当年 CSV 值；
# climate=:fair 时，T_global 已在 get_model 中连接到 temperature:T 组件输出，
# 此处不覆盖，T0_global 传 NaN → claim_damage 运行时取 damage_start_year
# 当年的 FaIR 温度作参考（SCC 计算时由 helper_scc.jl 显式覆写为 Base/Pulse
# 共用参考温度）。
function update_claim_damage_params!(m; data_dir::AbstractString, years, regions, theta_id::Int, start_year::Int, bau_ref_year::Int, climate::Symbol=:csv,
                                   damage_function::Symbol=:BK_region, income_state::Symbol=:dynamic, capital_response::Bool=true)
    shared_dir = joinpath(data_dir, "shared_parameters")
    unshared_dir = joinpath(data_dir, "unshared_parameters")
    centers = load_damage_state_centers(unshared_dir)
    basis = load_basis_cache(unshared_dir)
    theta_payload = load_theta_draws(unshared_dir)
    controls = claim_damage_controls(damage_function, income_state, capital_response)

    if climate == :csv
        T_global = _load_time_vector(joinpath(shared_dir, "climate_T_global.csv"), years)
        update_param!(m, :claim_damage, :T_global, T_global)
        update_param!(m, :claim_damage, :T0_global, T_global[findfirst(==(bau_ref_year), years)])
    else
        update_param!(m, :claim_damage, :T0_global, NaN)
    end
    update_param!(m, :claim_damage, :model_start_year, start_year)
    update_param!(m, :claim_damage, :damage_start_year, bau_ref_year)
    update_param!(m, :claim_damage, :state_window, BKW_STATE_WINDOW)
    # 局地渠道吃纯局地天气，全球年际变率单走全球渠道（open_issues.md P1-19）。
    # 两个文件都由 helper_mesmer.jl::regenerate_national_mesmer_inputs 与
    # climate_region_lv.csv 同一次抽样、同一 seed 生成，三者必须配套使用。
    update_param!(m, :claim_damage, :region_lv_local,
                  _load_time_region_matrix(joinpath(unshared_dir, "climate_region_lv_local.csv"), years, regions))
    update_param!(m, :claim_damage, :global_lv,
                  _load_time_vector(joinpath(unshared_dir, "climate_global_lv.csv"), years))
    update_param!(m, :claim_damage, :barT_center, centers["barT_center"])
    update_param!(m, :claim_damage, :lny0_center, centers["lny0_center"])
    update_param!(m, :claim_damage, :basis_Y1_g, matrix_from_channel(basis, regions, "G", "Y1"))
    update_param!(m, :claim_damage, :basis_Y2_g, matrix_from_channel(basis, regions, "G", "Y2"))
    update_param!(m, :claim_damage, :basis_Y1_nocap_g, matrix_from_channel(basis, regions, "G", "Y1_nocap"))
    update_param!(m, :claim_damage, :basis_Y2_nocap_g, matrix_from_channel(basis, regions, "G", "Y2_nocap"))
    update_param!(m, :claim_damage, :basis_Z1_g, matrix_from_channel(basis, regions, "G", "Z1"))
    update_param!(m, :claim_damage, :basis_Z2_g, matrix_from_channel(basis, regions, "G", "Z2"))
    update_param!(m, :claim_damage, :basis_Y1_l, matrix_from_channel(basis, regions, "L", "Y1"))
    update_param!(m, :claim_damage, :basis_Y2_l, matrix_from_channel(basis, regions, "L", "Y2"))
    update_param!(m, :claim_damage, :basis_Y1_nocap_l, matrix_from_channel(basis, regions, "L", "Y1_nocap"))
    update_param!(m, :claim_damage, :basis_Y2_nocap_l, matrix_from_channel(basis, regions, "L", "Y2_nocap"))
    update_param!(m, :claim_damage, :basis_Z1_l, matrix_from_channel(basis, regions, "L", "Z1"))
    update_param!(m, :claim_damage, :basis_Z2_l, matrix_from_channel(basis, regions, "L", "Z2"))

    models = Array{Float64}(undef, length(regions), 2)
    cmax = Array{Float64}(undef, length(regions), 2)
    for (i, r) in enumerate(regions)
        models[i, 1] = basis[string(r)]["G"]["model"]
        models[i, 2] = basis[string(r)]["L"]["model"]
        # C 的峰高界逐国不同（∝ D^ν），由 precompute_basis.py 预计算下发。
        cmax[i, 1] = basis[string(r)]["G"]["C_max"]
        cmax[i, 2] = basis[string(r)]["L"]["C_max"]
    end
    update_param!(m, :claim_damage, :basis_model, models)
    update_param!(m, :claim_damage, :basis_C_max, cmax)
    update_param!(m, :claim_damage, :w_g, weights_from_cache(basis, regions, "w_g"))
    update_param!(m, :claim_damage, :w_l, weights_from_cache(basis, regions, "w_l"))
    update_param!(m, :claim_damage, :temperature_effect, controls[:temperature_effect])
    update_param!(m, :claim_damage, :income_effect, controls[:income_effect])
    update_param!(m, :claim_damage, :capital_effect, controls[:capital_effect])

    for (name, value) in theta_parameters(theta_payload, theta_id, regions)
        update_param!(m, :claim_damage, name, value)
    end
    update_param!(m, :claim_damage, :alpha_c, [basis[string(r)]["G"]["alpha"] for r in regions])
    update_param!(m, :claim_damage, :delta_c, [basis[string(r)]["G"]["delta"] for r in regions])
    update_param!(m, :claim_damage, :rho_c, [basis[string(r)]["G"]["rho"] for r in regions])
end
