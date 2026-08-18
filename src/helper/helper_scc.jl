# ============================================================================
# 【helper_scc.jl】CLAIM 损害模块专用的 SCC 计算模块（取代已删除的 MimiFUND 移植文件
# new_marginaldamages.jl）。
# 设计要点：
#   1. 仅支持 CO₂ 脉冲，注入 FaIR 的 co2_cycle（模型须以 climate=:fair 构建）；
#   2. 年份索引全部基于模型自身的 time labels，不再硬编码 1990/1950 基年；
#   3. Base/Pulse 共用参考温度：先跑一次 Base，把 T0_global 的 NaN 哨兵
#      解析为 damage_start_year 当年的 FaIR 温度并显式覆写（coupling_plan §1.6），
#      避免脉冲污染 Pulse 轨的参考温度；
#   4. Ramsey 折现率 = prtp + eta·g，g 用含气候损害的 Base 轨
#      claim_damage:ypc_new 增长率逐国逐年递推（DICE 式内生折现——损害压低
#      未来增长 → 折现率内生下降 → SCC 上升）。
#   5. θ 不确定性不在此处抽样——外部循环 get_model(theta_id=i) 后逐次调用
#      compute_scco2 即可（damage_theta_draws.csv 预抽样机制）。
# ============================================================================

# 【emissionspulse 组件】排放脉冲加法器：output = input + add。
# 插在 emissions 组件与 co2_cycle 之间，用于给 Pulse 轨注入边际排放。
# 与 Mimi.adder 类似，但允许 missing 值直接透传。
@defcomp emissionspulse begin
    add    = Parameter(index=[time])
    input  = Parameter(index=[time])
    output = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output[t] = Mimi.@allow_missing(p.input[t]) + p.add[t]
    end
end

# 【_time_index】模型时间轴上 year 的下标；year 不在时间轴内时报错。
function _time_index(m, year::Int)
    idx = findfirst(==(year), Mimi.time_labels(m))
    idx === nothing && error("year $year 不在模型时间轴 $(first(Mimi.time_labels(m))):$(last(Mimi.time_labels(m))) 内")
    return idx
end

function active_damage_component(m::Model)
    Mimi.has_comp(m, :claim_damage) && return :claim_damage
    Mimi.has_comp(m, :bhm_damage) && return :bhm_damage
    error("model has no supported damage component")
end

# 【add_co2_pulse!】向模型 m 插入 emissionspulse 组件并接线：断开
# emissions 与 co2_cycle 的直连，改为经过脉冲加法器；把 pulse_size 吨 CO₂
# 平摊到 year 起的 pulse_length 年（换算为 GtC/yr：/1e9 吨→Gt，×12/44 CO₂→C）。
function add_co2_pulse!(m::Model, year::Int; pulse_size::Float64=1e7, pulse_length::Int=10)
    years = Mimi.time_labels(m)
    yi = _time_index(m, year)
    yi + pulse_length - 1 <= length(years) ||
        error("脉冲窗口 $(year)-$(year + pulse_length - 1) 超出模型时间轴末年 $(last(years))；请缩短 pulse_length 或延长 end_year")

    add_comp!(m, emissionspulse, after=:emissions)
    addem = zeros(length(years))
    addem[yi:(yi + pulse_length - 1)] .= pulse_size / (pulse_length * 1e9) * (12.0 / 44.0)
    update_param!(m, :emissionspulse, :add, addem)
    connect_param!(m, :emissionspulse, :input, :emissions, :E_co2)
    connect_param!(m, :co2_cycle, :E_co2, :emissionspulse, :output)
    return m
end

# 【get_marginal_model】构建 Base/Pulse 双轨边际模型对：以 m 为 Base，
# 复制出 Pulse 轨并注入 CO₂ 脉冲。两轨差分（除以 pulse_size）即得
# 单位排放的边际效应。
function get_marginal_model(m::Model; year::Int, pulse_size::Float64=1e7, pulse_length::Int=10)
    mm = create_marginal_model(m, pulse_size)
    add_co2_pulse!(mm.modified, year; pulse_size=pulse_size, pulse_length=pulse_length)
    return mm
end

# 【compute_scco2】计算指定排放年份的 CO₂ 社会成本（$/吨 CO₂，价格基准
# 与 ypc 数据一致，即 PWT 2017 不变价美元）。
# 流程：(1) 跑 Base 解析共用参考温度 T0；(2) 构建 Base/Pulse 边际模型对并
# 运行至 last_year；(3) 差分得边际损害（claim_damage:loss，$/吨），按
# Ramsey 折现（prtp + eta·g）逐国折现后加总。
# 关键字参数：
#   year         —— 排放脉冲年份（必填）
#   last_year    —— 损害积分截断年，默认取模型时间轴末年（get_model 默认
#                   2150；时间轴最远可建到 2500，用户可任选截断年）
#   eta, prtp    —— Ramsey 折现参数（默认 1.45 / 0.015，FUND 口径）
#   pulse_size   —— 脉冲总量（吨 CO₂，默认 1e7）
#   pulse_length —— 脉冲平摊年数（默认 10，FUND 惯例）
#   return_details —— true 时返回 (scc, scc_regional, marginaldamage, df, mm)
function compute_scco2(m::Model; year::Int, last_year::Union{Int,Nothing}=nothing,
                       eta::Float64=1.45, prtp::Float64=0.015,
                       pulse_size::Float64=1e7, pulse_length::Int=10,
                       return_details::Bool=false)
    Mimi.has_comp(m, :co2_cycle) ||
        error("compute_scco2 需要 FaIR 气候链：请用 get_model(:national; climate=:fair) 构建模型")

    years = Mimi.time_labels(m)
    last_year = last_year === nothing ? Int(last(years)) : last_year
    li = _time_index(m, last_year)
    yi = _time_index(m, year)
    yi <= li || error("排放年份 year=$year 须不晚于 last_year=$last_year")

    # (1) Base 先行运行，解析 Base/Pulse 共用的 B&K 参考温度 T0
    run(m)
    damage_component = active_damage_component(m)
    if damage_component == :claim_damage
        T0 = m[:claim_damage, :T0_global]
        if isnan(T0)
            dsy = Int(m[:claim_damage, :damage_start_year])
            T0 = m[:temperature, :T][_time_index(m, dsy)]
            update_param!(m, :claim_damage, :T0_global, T0)
        end
    end

    # (2) 构建并运行边际模型对（Base 以显式参考温度重跑，Pulse 含排放脉冲）
    mm = get_marginal_model(m; year=year, pulse_size=pulse_size, pulse_length=pulse_length)
    run(mm; ntimesteps=li)

    # (3) 边际损害与 Ramsey 折现。population_sspX.csv 为 PWT 口径（百万人），
    # 故 loss 单位是百万美元，差分后 ×1e6 换算为 $/吨 CO₂。
    marginaldamage = mm[damage_component, :loss] .* 1e6
    ypc = mm.base[damage_component, :ypc_new]
    nregions = size(marginaldamage, 2)

    # 期次约定与 MimiFUND 一致（new_marginaldamages.jl:246-254）：**先赋值再前推**，
    #   df_τ = 1,  df_{t+1} = df_t / (1 + prtp + eta·g_t),  g_t = ypc_t/ypc_{t-1} − 1
    # 即 t 年的折现因子只吃到 t−1 年为止的增长。此前写成
    # df_t = df_{t-1}/(1 + prtp + eta·g_t)，整条折现路径比 FUND 晚一期；
    # 比值为 (1+prtp+eta·g_τ)/(1+prtp+eta·g_t)，g 随时间下降时会系统性少折现。
    #
    # ypc 取 mm.base 的 ypc_new——**不打脉冲**那一轨的**含损害**人均 GDP。
    # 两条轨必须共用同一条折现因子：各用各的，(5.11) 的差分里就会混进
    # "脉冲把折现率推动了多少"，那不是边际损害。
    yi > 1 || error("脉冲年份须晚于模型首年：折现率要用**进入**该年的增长率 ypc_t/ypc_{t-1}")
    df = zeros(li, nregions)
    for r in 1:nregions
        x = 1.0
        for t in yi:li
            df[t, r] = x
            g = ypc[t, r] / ypc[t - 1, r] - 1.0
            isfinite(g) || (g = 0.0)
            denom = 1.0 + prtp + eta * g
            denom > 0.0 || error("折现因子分母非正（region=$r, t=$t, g=$g）：请检查情景数据")
            x = x / denom
        end
    end

    scc_regional = vec(sum(marginaldamage[yi:li, :] .* df[yi:li, :], dims=1))
    scc = sum(scc_regional)

    if return_details
        return (scc=scc, scc_regional=scc_regional, marginaldamage=marginaldamage, df=df, mm=mm)
    end
    return scc
end
