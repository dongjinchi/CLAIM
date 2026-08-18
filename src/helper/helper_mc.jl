# ============================================================================
# 【helper_mc.jl】外层气候不确定性 Monte Carlo（open_issues.md P1-7）
#
# 气候模块的不确定性分两段，本模块同时抽样：
#   ① 全球平均温升响应 —— FaIR 参数集合 climate_fair_ensemble.csv 的成员
#      （热响应 q/d、CO₂ 强迫 f、碳循环 r、气溶胶 ari/aci），由
#      helper_fairv2.jl::apply_fair_member! 应用。该函数与构建集合时的打分
#      共用同一实现，见其文件内注释。**按 seed 随机置换后取前 n 行**，理由见
#      run_climate_mc 内 fair_idx 处的注释——集合行序有含义，取前缀是有偏的。
#   ② 降尺度的**升温模式** —— MESMER v1.0.0 的 58 个 CMIP6 模式，逐抽样换 ESM，
#      由 apply_region_lt! 覆写国家级 β。
#   ③ 降尺度的内部变率 —— MESMER 局地变率实现，逐抽样换 seed。
#
# ②③ 是 beusch2020mesmer 明确区分的两层：模式间散布在年尺度上 "clearly exceeds"
# 单模式初始条件散布。本项目实测印证了这个次序——③ 只值 SCC 的 0.19%，而 ② 的
# 逐国 β 变异系数中位就有 10.8%（open_issues.md 结构性 caveat #7）。
#
# 关于 ③ 的实现方式：**每次抽样现场在网格上模拟**，不从预生成的实现池里取。
# 池化会退回 P1-3 之前「读现成 emu 文件」的模式，使 helper_mesmer.jl 的网格
# 模拟形同虚设；本项目的目标是模型自身端到端跑通 FaIR 与 MESMER，而非外接
# 参数。实测代价也不支持池化：网格缓存后约 2.6 s/次抽样。
#
# 损害参数 θ 的抽样是**内层**，不在本模块：外部循环 get_model(theta_id=i) 即可
# （damage_theta_draws.csv 预抽样机制，指南 §3.5.1）。本模块固定 theta_id，
# 以便把气候维度的贡献单独分离出来。
# ============================================================================

using CSV, DataFrames, Mimi, Printf, Random, Statistics

# 【apply_region_lv!】把一条局地变率实现写进已建好的模型（内存内覆写）。
#
# 与既有加载路径的分工——两者职责不重叠，都要保留：
#   * 构造时：get_model 经 helper_global.jl / helper_damage.jl 从磁盘的
#     climate_region_lv.csv 读入。生产运行与单跑走这条路。
#   * 运行时：本函数用内存中现场模拟的矩阵覆写。仅 MC 与反事实实验走这条路，
#     避免每次抽样都重写 CSV。
#
# 两个必须守住的不变量：
#   1. **必须在 compute_scco2 之前调用**。SCC 的 Base/Pulse 两轨要共用同一条 ν，
#      否则局地渠道不再在边际中相消（指南 §3.6）。
#   2. **消费 region_lv 的组件全部要覆写**：climateregional 用它算 ts→regtmp 再
#      喂给损害组件，而损害组件又单独用它走局地渠道。只写其中一处会让两者不
#      自洽，而且**不会报错**。故此处用 active_damage_component 判断当前活跃的
#      损害组件，不硬编码 :claim_damage——那样会静默跳过 :bhm_damage。
# 第三个不变量（P1-19）：climateregional 拿完整 ν（决定国别温度产出），而
# 损害组件的局地渠道只能拿**纯局地**分量，全球年际变率单走 global_lv。
# 三者来自同一次抽样，故本函数一律接 simulate_region_lv_components 的返回值，
# 不接单个矩阵——分开传极易把完整 ν 误接到局地渠道上。
function apply_region_lv!(m::Model, c::NamedTuple)
    update_param!(m, :climateregional, :region_lv, c.total)
    dmg = active_damage_component(m)
    if dmg === :claim_damage
        update_param!(m, dmg, :region_lv_local, c.loc)
        update_param!(m, dmg, :global_lv, c.gv)
    else
        # bhm_damage 用的是绝对温度基线 tempbase + region_lv，不做冲击分解，
        # 因此仍消费完整 ν。
        update_param!(m, dmg, :region_lv, c.total)
    end
    return nothing
end

# 【load_lt_ensemble】加载 58 个 CMIP6 模式的国家级升温模式系数。
#
# 数据来自 build_mesmer_v1_ensemble.py（MESMER v1.0.0 预标定参数 W·β 投影），
# 与自训练的 climate_ltcoef.csv 是**两套网格**：集合走 g025 2957 格点，自训练
# 走 IPSL 原生 5940 格点。两者不可按下标混用，故各自带自己的权重矩阵。
#
# 返回矩阵按 regions 的顺序重排，不依赖 CSV 的行序；任何一个国家缺失都会因
# NaN 检查而报错，不静默补零——旧链路的静默零值正是 P1-3 的病灶。
function load_lt_ensemble(data_dir::AbstractString, regions;
                          path::AbstractString="")
    p = isempty(path) ?
        joinpath(data_dir, "unshared_parameters", "climate_ltcoef_ensemble.csv") : path
    df = CSV.read(p, DataFrame)
    esms = unique(String.(df.esm))
    eidx = Dict(e => i for (i, e) in enumerate(esms))
    ridx = Dict(String(r) => i for (i, r) in enumerate(regions))
    coef = fill(NaN, length(esms), length(regions))
    inter = fill(NaN, length(esms), length(regions))
    for row in eachrow(df)
        r = get(ridx, String(row.region), 0)
        r == 0 && continue
        coef[eidx[String(row.esm)], r] = row.ltcoef
        inter[eidx[String(row.esm)], r] = row.ltinter
    end
    if any(isnan, coef) || any(isnan, inter)
        i, j = Tuple(findfirst(isnan, coef))
        error("climate_ltcoef_ensemble.csv 缺 $(esms[i]) × $(regions[j])")
    end
    return (esms=esms, coef=coef, inter=inter,
            weight=_load_esm_weights(dirname(p), esms))
end

# 【_load_esm_weights】读模式独立性权重（Knutti et al. 2017），与集合同目录。
# 文件缺失时退回等权并给出提示——等权是可用的对照配置，不是错误。
function _load_esm_weights(dir::AbstractString, esms)
    path = joinpath(dir, "climate_ltcoef_ensemble_weights.csv")
    isfile(path) || return fill(1 / length(esms), length(esms))
    df = CSV.read(path, DataFrame)
    lookup = Dict(String(r.esm) => Float64(r.weight) for r in eachrow(df))
    w = [get(lookup, e, NaN) for e in esms]
    any(isnan, w) && error("climate_ltcoef_ensemble_weights.csv 缺模式 " *
                           join(esms[isnan.(w)], ", "))
    all(w .> 0) || error("模式权重必须为正")
    return w ./ sum(w)
end

# 【_sample_categorical】按给定概率抽 n 个下标。逆变换法，不引入 StatsBase。
function _sample_categorical(rng, probs::AbstractVector, n::Int)
    cum = cumsum(probs)
    cum[end] = 1.0                       # 防浮点累积误差让 rand() 落在末尾之外
    return [searchsortedfirst(cum, rand(rng)) for _ in 1:n]
end

# 【apply_region_lt!】把某个 ESM 的国家级升温模式写进已建好的模型。
#
# 与 apply_region_lv! 的关键差别：**region_lt_coef 只有 climateregional 一个
# 消费者**（已核查 src 全域），故不需要 active_damage_component 那套分发；
# 损害组件拿到的是 climateregional 算好的 regtmp，不直接读 β。
#
# region_lt_inter（γ）在生产配置下**不进入前向段**：水平锚由 30 年观测窗口给出，
# 截距只在无观测数据的 legacy 分支与诊断里用（ClimateRegionalComponent.jl 顶部
# 注释）。这里仍一并覆写，是为了让两个分支自洽——只写其一不会报错，但会让
# legacy 分支拿到 A 模式的斜率配 B 模式的截距。
function apply_region_lt!(m::Model, coef::AbstractVector, inter::AbstractVector)
    update_param!(m, :climateregional, :region_lt_coef, collect(Float64, coef))
    update_param!(m, :climateregional, :region_lt_inter, collect(Float64, inter))
    return nothing
end

# 【draw_region_lv】为第 seed 次抽样现场模拟一条国家级变率实现。
#
# gen_years 必须与生产文件 climate_region_lv.csv 的生成区间一致（1950:2500），
# 之后再截到模型年限。原因：RNG 的抽取次数依赖模拟长度，若只在模型年限上模拟，
# **同一 seed 会得到不同的实现**，与生产基准对不上。截断后 seed=42 与生产
# 逐值相同。
function draw_region_lv(grid, model_years, seed::Integer;
                        gen_years=collect(1950:2500))
    c = simulate_region_lv_components(grid, gen_years; seed=seed)
    n = length(model_years)
    return (total=c.total[1:n, :], loc=c.loc[1:n, :], gv=c.gv[1:n])
end

# 【apply_theta!】把第 theta_id 次 pooled θ 抽样写进已建好的模型（内存内覆写）。
#
# 与 get_model(theta_id=i) 的分工，同 apply_region_lv! 那一组：构造时走
# update_claim_damage_params!，运行时走本函数。**这正是联合抽样能做的原因**——
# get_model 要 2.84 s 而 run 只要 0.17 s，逐抽样重建模型会把 20 分钟变成 3 小时。
#
# θ 只有 6 个 pooled 系数（各展开成 n_regions × 11 的矩阵），国家异质性由运行时
# 的 T_state/Y_state 投影内生产生，故此处不需要按国分发。
#
# 只对 :claim_damage 有意义：bhm_damage 用的是 BHM 二次型，没有 θ 这一层。抽了
# 也不会报错、只是完全无效，故显式拒绝而不是静默跳过。
function apply_theta!(m::Model, theta_payload, theta_id::Int, regions)
    dmg = active_damage_component(m)
    dmg === :claim_damage ||
        error("apply_theta! 只适用于 :claim_damage，当前活跃损害组件是 $(dmg)")
    for (name, value) in theta_parameters(theta_payload, theta_id, regions)
        update_param!(m, :claim_damage, name, value)
    end
    return nothing
end

# 【_mc_mode_label】由四个开关**导出**本次运行的口径标签，而不是另设一个 mode 参数。
#
# 设成参数就会出现"mode 写着 :joint、开关却只开了气候"这种自相矛盾且不报错的状态。
# P1-12 记录过三次同型教训：声明的意图与实现只在特例下等价。标签是对实际抽了什么
# 的描述，不是控制项。
function _mc_mode_label(sample_fair::Bool, sample_esm::Bool, sample_mesmer::Bool,
                        sample_theta::Bool)
    climate = sample_fair || sample_esm || sample_mesmer
    climate && sample_theta && return "joint"
    climate && return "climate"
    sample_theta && return "damage"
    return "none"
end

# 【run_scc_mc】SCC 的外层不确定性 MC 主循环（气候 × 损害）。
#
# 2026-08-06 由 run_climate_mc 推广而来：θ 从"外部循环重建模型"改为第四个开关
# sample_theta，与既有三个开关同构（open_issues.md §6.3 H2 / R24）。
#
# **四个维度各用独立 RNG，这不是风格问题而是设计要求**：只有这样，关掉任一维时
# 其余维的抽样序列才原样不动，三个口径（气候边际 / 损害边际 / 联合）才共用同一
# 串随机数。这就是共同随机数（common random numbers）——joint 与 climate 之差
# 于是干净地只反映 θ 层，抽样噪声大部分相消，比独立跑三次精确得多。
#
# 逐次追加写出（append=k>1），因此任何时候中断都保得住已完成的部分，也支持
# 边跑边看收敛（coupling_plan.md §9 的「MC 收敛检验」正需要这种形态）。
#
# 模型只建一次，之后原地 update_param! 换参数：实测 get_model 要 2.84 s 而
# run 仅 0.17 s，反复重建会把耗时放大一个量级。
function run_scc_mc(; n::Int=200,
                        scenario::String="ssp245",
                        pulse_year::Int=2030,
                        last_year::Int=2100,
                        end_year::Int=2150,
                        theta_id::Int=0,
                        data_dir::AbstractString=joinpath(@__DIR__, "..", "..", "data"),
                        ensemble_path::AbstractString="",
                        output_path::AbstractString="",
                        sample_fair::Bool=true,
                        sample_esm::Bool=true,
                        sample_mesmer::Bool=true,
                        lt_ensemble_path::AbstractString="",
                        esm_weighting::Symbol=:independence,
                        esm_seed::Int=20260803,
                        fair_sampling::Symbol=:random,
                        fair_seed::Int=20260805,
                        sample_theta::Bool=false,
                        theta_seed::Int=20260806,
                        verbose::Bool=true)
    data_dir = abspath(data_dir)
    ens_path = isempty(ensemble_path) ?
               joinpath(data_dir, "unshared_parameters", "climate_fair_ensemble.csv") :
               ensemble_path
    mode = _mc_mode_label(sample_fair, sample_esm, sample_mesmer, sample_theta)
    # 默认输出路径按口径分开。三个口径写同一个文件会互相覆盖，而且覆盖不报错——
    # 默认值本身就应当挡住这件事，不能指望调用方每次记得传 output_path。
    out_path = isempty(output_path) ?
               abspath(joinpath(@__DIR__, "..", "..", "results", "scc_mc_$(mode).csv")) :
               output_path

    ens = CSV.read(ens_path, DataFrame)
    n <= nrow(ens) || error("n=$n exceeds ensemble size $(nrow(ens))")

    # FaIR 成员：**必须随机置换后取，不能取前 n 行**（open_issues.md §6.2 R7 / P1-7）。
    #
    # climate_fair_ensemble.csv 是 build_fair_constrained_ensemble.py 的
    # systematic_resample 产物：2000 行等权，权重由**重复次数**编码，行序即源先验的
    # 生成顺序（source_member_id 单调不减）。因此任何连续切片都是聚簇样本，不是子样本：
    # 前 300 行只覆盖 **119 个不同**的源成员，而随机 300 行平均覆盖 241 个
    # （200 次重抽的范围 226–255）——有效样本量少一半。分布上也可区分，前 300 行
    # 对其余 1700 行的两样本 KS：ECS p=2.3e-4、TCR p=1.7e-5、rT p=1.2e-7，且三者
    # 均偏暖（ECS 2.710 vs 2.658、TCR 1.792 vs 1.754、rT 2.499 vs 2.246）。
    #
    # **按行等概率抽，不要去重**：重复正是后验权重的编码方式，去重等于把加权后验
    # 拍回等权先验。
    #
    # 独立 RNG，理由同下面的 esm_rng。:prefix 保留旧行为，仅用于复现 2026-08-05
    # 之前的产物，不可用于生产。
    fair_idx = if fair_sampling === :random
        randperm(MersenneTwister(fair_seed), nrow(ens))[1:n]
    elseif fair_sampling === :prefix
        collect(1:n)
    else
        error("fair_sampling 只能是 :random 或 :prefix，收到 $(fair_sampling)")
    end
    if verbose && sample_fair
        @printf("FaIR 抽样: %s，%d 行集合中取 %d 行（覆盖 %d 个不同源成员）\n",
                fair_sampling, nrow(ens), n,
                length(unique(ens.source_member_id[fair_idx])))
    end

    regions = load_national_regions(joinpath(data_dir, "shared_parameters"))

    # 损害 θ：读 damage_theta_draws.csv 的预抽样，按 draw_id 取，Julia 侧不生成
    # 随机数（指南 §3.5.1）。与 ①②③ 不同的理由记在 CLAIM_model_guide §4.7：
    # 抽样在 Python 侧用 numpy PCG64 完成，预生成才能跨语言逐位复现，也才能被
    # test_theta_draws_stats.py 逐条审计。
    #
    # **与 FaIR 集合的关键差别**：这里每行是一次独立的 MVN 抽样，**没有用重复
    # 次数编码权重**，所以不放回抽样只是效率略优，不涉及后验失真。不要把 R7 的
    # 推理照搬过来。
    theta_payload = sample_theta ?
                    load_theta_draws(joinpath(data_dir, "unshared_parameters")) : nothing
    theta_idx = if sample_theta
        n_theta = length(theta_payload["theta_draws"])
        n_theta > 0 || error("damage_theta_draws.csv 为空，无法抽样 θ")
        n <= n_theta ||
            error("n=$(n) 超过 θ 抽样数 $(n_theta)，请用 precompute_mvn_draws.py --n_draws 重新生成")
        randperm(MersenneTwister(theta_seed), n_theta)[1:n]
    else
        fill(theta_id, n)
    end
    verbose && sample_theta &&
        @printf("θ 抽样: %d 条预抽样中取 %d 条（seed %d）\n",
                length(theta_payload["theta_draws"]), n, theta_seed)

    m = get_model(:national; climate=:fair, scenario=scenario, end_year=end_year,
                  theta_id=theta_id, data_dir=data_dir)
    years = collect(Mimi.time_labels(m))
    grid = sample_mesmer ? cached_mesmer_grid(data_dir, regions) : nothing

    # ESM 用独立 RNG 抽，不复用 draw 序号：这样开关 sample_esm 不会改变 FaIR
    # 成员与变率 seed 的配对，三个维度可以各自单独开关做归因。
    #
    # 默认按**模式独立性权重**抽（Knutti et al. 2017，标定见
    # build_mesmer_v1_ensemble.py）。等权抽 58 个模式会放大 CMIP6 的模式家族
    # ——EC-Earth3 五个变体、E3SM 四个、CanESM5/CESM2/GISS/CNRM 各三个——
    # 等于让同一个模式投多票。esm_weighting=:equal 保留等权作为对照配置。
    lt = sample_esm ? load_lt_ensemble(data_dir, regions; path=lt_ensemble_path) : nothing
    esm_rng = MersenneTwister(esm_seed)
    esm_draws = if !sample_esm
        fill(0, n)
    elseif esm_weighting === :independence
        _sample_categorical(esm_rng, lt.weight, n)
    elseif esm_weighting === :equal
        rand(esm_rng, 1:length(lt.esms), n)
    else
        error("esm_weighting 只能是 :independence 或 :equal，收到 $(esm_weighting)")
    end
    if verbose && sample_esm
        @printf("ESM 抽样: %s 加权，%d 个模式中抽到 %d 个（有效模式数 %.1f）\n",
                esm_weighting, length(lt.esms), length(unique(esm_draws)),
                1 / sum(abs2, lt.weight))
    end

    base_scc = Float64(compute_scco2(m; year=pulse_year, last_year=last_year))
    verbose && @printf("基准 SCC(%d, %s, →%d) = %.2f \$/tCO2\n\n",
                       pulse_year, scenario, last_year, base_scc)

    mkpath(dirname(out_path))
    # source_member_id 单列出来，是因为它才是"抽到了几个不同的气候参数集"的口径：
    # member_id 是重采样后的行号，逐行互异，看不出重复。R7 正是靠这一列发现的。
    # mode 与 theta_id 落盘，三个口径的产物才能直接纵向拼接做方差分解。
    # 未抽样的维度一律写 -1 / NaN 而不是写"当时恰好是什么"：损害边际里若照抄
    # 一个没被应用的 FaIR 成员的 ECS/TCR，读表的人会以为那次抽样用了它。
    rows = DataFrame(draw=Int[], mode=String[], member_id=Int[], source_member_id=Int[],
                     esm=String[], seed=Int[], theta_id=Int[],
                     ECS=Float64[], TCR=Float64[], beta_med=Float64[],
                     scc=Float64[], T2100=Float64[], secs=Float64[])
    for k in 1:n
        t0 = time()
        row = ens[fair_idx[k], :]
        sample_fair && apply_fair_member!(m, row)
        e = esm_draws[k]
        sample_esm && apply_region_lt!(m, view(lt.coef, e, :), view(lt.inter, e, :))
        sample_mesmer && apply_region_lv!(m, draw_region_lv(grid, years, k))
        sample_theta && apply_theta!(m, theta_payload, theta_idx[k], regions)
        run(m)
        T2100 = Float64(m[:temperature, :T][findfirst(==(2100), years)])
        scc = Float64(compute_scco2(m; year=pulse_year, last_year=last_year))
        rec = (k, mode,
               sample_fair ? Int(row.member_id) : -1,
               sample_fair ? Int(row.source_member_id) : -1,
               sample_esm ? lt.esms[e] : "baseline",
               sample_mesmer ? k : -1,
               theta_idx[k],
               sample_fair ? Float64(row.ECS) : NaN,
               sample_fair ? Float64(row.TCR) : NaN,
               sample_esm ? median(view(lt.coef, e, :)) : NaN,
               scc, T2100, time() - t0)
        push!(rows, rec)
        # 逐次落盘。写 rows 的最后一行而非 DataFrame([rec])：后者从裸 Tuple 建表
        # 会把列名退化成 1..10，落盘文件因此丢失全部列名且不报错。
        CSV.write(out_path, rows[end:end, :]; append=(k > 1))
        # 按列名取，不按元组下标：加一列就会让下标全部错位，而 %s 打印一个 Int
        # 不会报错，只会静默印错东西。
        verbose && (k <= 3 || k % 50 == 0) &&
            @printf("  draw %4d/%d  ECS %5.2f  TCR %5.2f  %-18s θ%-5d SCC %8.2f  (%.1f s)\n",
                    k, n, rows.ECS[end], rows.TCR[end], rows.esm[end],
                    rows.theta_id[end], scc, rows.secs[end])
    end

    if verbose
        v = sort(rows.scc)
        q(p) = v[clamp(ceil(Int, p * length(v)), 1, length(v))]
        @printf("\n[%s] %d 次抽样  用时合计 %.1f 分钟（中位 %.1f s/次）\n",
                mode, n, sum(rows.secs) / 60, median(rows.secs))
        @printf("SCC  5%% %.2f   50%% %.2f   95%% %.2f   均值 %.2f   sd %.2f\n",
                q(0.05), q(0.50), q(0.95), mean(v), std(v))
        @printf("相对基准 %.2f: [%+.1f%%, %+.1f%%]   MCSE(均值) %.2f\n",
                base_scc, 100 * (q(0.05) / base_scc - 1), 100 * (q(0.95) / base_scc - 1),
                std(v) / sqrt(length(v)))
        println("写入 ", out_path)
    end
    return rows
end

# 【run_climate_mc】历史入口名，保留以免旧脚本与文档片段失效。
# 语义等价于 run_scc_mc(; sample_theta=false, ...)，即"只抽气候"这一个口径。
# 新代码请直接用 run_scc_mc，它的名字才反映实际抽了什么。
run_climate_mc(; kwargs...) = run_scc_mc(; sample_theta=false, kwargs...)
