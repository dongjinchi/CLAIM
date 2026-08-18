# ============================================================================
# 【claim_damage 组件】B&K（Bilal & Känzig）损害核 Mimi 组件——本项目的新版损害模块。
# 功能：每年每国执行"动态θ投影 → 闭式OLS反演(A,C) → 组装损害核ζ → 温度卷积得TFP
#       损害 → Solow资本递归 → 映射回真实人均GDP并计算损失"的完整链条。
# 理论依据：coupling_plan.md §1.4-1.6；CLAIM_model_guide.md §3.4。
# 关键设计：TFP 卷积仅 11 年记忆（ζ 指数衰减），资本 khat 为状态变量、无限记忆；
#           全球渠道用温度水平偏离(T_global - T0_global)，局部渠道直接用平稳的
#           region_lv 天气噪声（不减参考年）。
#           T0_global 传 NaN 时在运行时取 damage_start_year 当年的 T_global
#           （FaIR 接入模式下建模时无法预知参考温度）；SCC 计算时应显式传入
#           Base/Pulse 共用的参考温度（coupling_plan §1.6）。
#
# 温度状态的时间口径（open_issues.md P1-4，2026-07-30 修改）
# ---------------------------------------------------------------------------
# LP 中的 barT_c 是国家 1950–1979 的**三十年气候态均值**。此前 T_state 用的是
# 滞后一年的 regtmp 并扣掉 region_lv，"一年值"与"三十年均值"语义不同：接入
# FaIR + MESMER 后 ts_lt 自身的年际波动会直接进入 θ 的状态投影。现改为
# **过去 30 年 regtmp 的滚动均值**，窗口 [t-30, t-1]：
#
#   * 窗口长度 30 不是自由参数，它由 LP 端 barT_s 的 1950–1979 窗口钉死；
#   * **不扣 region_lv**。LP 的 barT_s 本身就是从实测 lctmp_bkly_pw 算的
#     三十年均值（含天气），WMO 的"气候常态"定义亦然，故含 ν 才是同源口径。
#     30 年平均已把噪声压到国家级 sd 的约 1/4（残留中位 0.077 °C，约为 barT
#     跨国 sd 7.49 °C 的 1%），不必再显式剔除。
#     代价：θ 因此依赖 MESMER 的随机种子，P1-7 外层气候 MC 的"气候 vs 损害
#     参数"分解会被轻微污染——这是刻意接受的取舍。
#   * SCC 不受影响：Base/Pulse 共用同一条 ν，滚动均值之差只剩受迫分量。
#
# θ 随气候态时变（而非冻结在 LP 的历史气候态）是刻意选择：θ₁ 若真是纯混淆
# 就不该在横截面上使用（那等于退回 damage_function=:BK_global），既然用了它
# 解释国家间差异，就必须承认同一梯度在时间上也生效。这属可迁移性假设，
# 见 open_issues.md §3 第 1 条。
# ============================================================================
@defcomp claim_damage begin
    regions = Index()
    horizon_dim = Index()

    T_global = Parameter(index=[time])
    T0_global = Parameter()
    regtmp = Parameter(index=[time, regions])
    # 局地渠道只吃**纯特异性局地天气**，不含全球变率驱动的分量；后者随
    # global_lv 并入全球渠道。两个参数的口径见下方 z_g / z_l 处的注释与
    # open_issues.md P1-19。刻意不叫 region_lv：climateregional 的 region_lv
    # 是完整 ν，同名不同义会静默算错，改名可让漏改的接线直接报错。
    region_lv_local = Parameter(index=[time, regions])
    global_lv = Parameter(index=[time])

    barT_center = Parameter()
    lny0_center = Parameter()

    ypc_baseline = Parameter(index=[time, regions])
    population = Parameter(index=[time, regions])

    basis_model = Parameter(index=[regions, 2])
    # C 的峰高界 0.37·e^ν·(D/ν)^ν，逐国逐渠道不同（open_issues.md P1-12）。
    basis_C_max = Parameter(index=[regions, 2])
    basis_Y1_g = Parameter(index=[regions, horizon_dim])
    basis_Y2_g = Parameter(index=[regions, horizon_dim])
    basis_Y1_nocap_g = Parameter(index=[regions, horizon_dim])
    basis_Y2_nocap_g = Parameter(index=[regions, horizon_dim])
    basis_Z1_g = Parameter(index=[regions, horizon_dim])
    basis_Z2_g = Parameter(index=[regions, horizon_dim])
    basis_Y1_l = Parameter(index=[regions, horizon_dim])
    basis_Y2_l = Parameter(index=[regions, horizon_dim])
    basis_Y1_nocap_l = Parameter(index=[regions, horizon_dim])
    basis_Y2_nocap_l = Parameter(index=[regions, horizon_dim])
    basis_Z1_l = Parameter(index=[regions, horizon_dim])
    basis_Z2_l = Parameter(index=[regions, horizon_dim])

    theta_g0 = Parameter(index=[regions, horizon_dim])
    theta_g1 = Parameter(index=[regions, horizon_dim])
    theta_g2 = Parameter(index=[regions, horizon_dim])
    theta_l0 = Parameter(index=[regions, horizon_dim])
    theta_l1 = Parameter(index=[regions, horizon_dim])
    theta_l2 = Parameter(index=[regions, horizon_dim])

    w_g = Parameter(index=[regions, horizon_dim])
    w_l = Parameter(index=[regions, horizon_dim])

    alpha_c = Parameter(index=[regions])
    delta_c = Parameter(index=[regions])
    # rho_c 是**传入的**资本留存权重，由 03_invert_damage.py 从实测资本存量
    # （PWT rnna）算出：rho_t = (1-delta)K_{t-1}/K_t，全样本均值。
    # 2026-08-08 之前这里声明的是 g_c 并现算 rho=(1-delta)/(1+g)——那是拿人均
    # GDP 增长率代理总资本增长率，实测单向偏长（半衰期 11.6 vs 7.9 年），
    # 且给了反演端与运行时两个真源。改名而不是改含义：漏改的接线会直接报错。
    rho_c = Parameter(index=[regions])
    temperature_effect = Parameter()
    income_effect = Parameter()
    capital_effect = Parameter()
    model_start_year = Parameter()
    damage_start_year = Parameter()
    # 温度状态的回看窗口长度（年）。由 LP 端 barT_s 的 1950–1979 三十年窗口
    # 钉死，不是自由参数；作为 Parameter 暴露只为让测试能覆盖短窗口边界。
    state_window = Parameter()

    khat = Variable(index=[time, regions])
    ypc_new = Variable(index=[time, regions])
    loss = Variable(index=[time, regions])
    yhat = Variable(index=[time, regions])
    z_tfp = Variable(index=[time, regions])
    clipped = Variable(index=[time, regions])

    # 【run_timestep】逐时间步执行损害计算：damage_start_year 之前输出基准值；
    # 之后按七步流程计算每国的 TFP 损害、资本偏离与实际损失（详见文件头注释）。
    function run_timestep(p, v, d, t)
        start_idx = max(1, Int(p.damage_start_year) - Int(p.model_start_year) + 1)
        # T0_global=NaN → 运行时取 damage_start_year 当年温度作参考（此时
        # temperature 组件已算完当期，读取安全；启动前的时间步不会进入该分支）。
        T0 = (t.t >= start_idx && isnan(p.T0_global)) ? Float64(p.T_global[TimestepIndex(start_idx)]) : p.T0_global

        # 四个长度 11 的工作数组提到 for r 之外复用（open_issues.md P2-9）。
        # 原先在循环内用推导式新建，185 国 × 201 年 × 4 个数组，实测占 run(m)
        # 全部分配的 68.5%（56.5 / 82.5 MB）。每个国家进入循环体后会先被完整
        # 覆写再使用，故跨国复用是安全的；提到 run_timestep 层而不用模块级全局，
        # 是为了不引入可变共享状态、也不依赖"Mimi 一定单线程"这个假设。
        nh = length(d.horizon_dim)
        th_g = Vector{Float64}(undef, nh)
        th_l = Vector{Float64}(undef, nh)
        zeta_g = Vector{Float64}(undef, nh)
        zeta_l = Vector{Float64}(undef, nh)

        for r in d.regions
            if t.t < start_idx
                v.khat[t, r] = 0.0
                v.ypc_new[t, r] = p.ypc_baseline[t, r]
                v.loss[t, r] = 0.0
                v.yhat[t, r] = 0.0
                v.z_tfp[t, r] = 0.0
                v.clipped[t, r] = 0.0
            else
                lag_idx = TimestepIndex(max(1, t.t - 1))
                # 气候态 = 过去 state_window 年 regtmp 的均值（窗口 [t-w, t-1]，
                # 不含当年，与"滞后状态"一致）。历史不足 w 年时用扩张窗口，
                # 绝不静默截断；生产设定下损害 2025 起算、模型 1950 起步，
                # 窗口始终是完整的 30 年。
                win_lo = max(1, t.t - Int(p.state_window))
                win_hi = max(1, t.t - 1)
                climatology = 0.0
                for i in win_lo:win_hi
                    climatology += p.regtmp[TimestepIndex(i), r]
                end
                climatology /= (win_hi - win_lo + 1)
                T_state = p.temperature_effect * (climatology - p.barT_center)
                lag_income = t.t == 1 ? p.ypc_baseline[lag_idx, r] : v.ypc_new[lag_idx, r]
                Y_state = p.income_effect * (log(lag_income) - p.lny0_center)

                # 显式计数器填缓冲区，逐元素与原推导式等价（推导式按迭代顺序
                # 产生元素，这里的 i 就是那个顺序）。
                i = 0
                for h in d.horizon_dim
                    i += 1
                    th_g[i] = p.theta_g0[r, h] + p.theta_g1[r, h] * T_state + p.theta_g2[r, h] * Y_state
                    th_l[i] = p.theta_l0[r, h] + p.theta_l1[r, h] * T_state + p.theta_l2[r, h] * Y_state
                end

                # 只读切片一律取视图，避免每国每年再复制 6 个长度 11 的数组。
                Y1_g = p.capital_effect > 0.5 ? (@view p.basis_Y1_g[r, :]) : (@view p.basis_Y1_nocap_g[r, :])
                Y2_g = p.capital_effect > 0.5 ? (@view p.basis_Y2_g[r, :]) : (@view p.basis_Y2_nocap_g[r, :])
                Y1_l = p.capital_effect > 0.5 ? (@view p.basis_Y1_l[r, :]) : (@view p.basis_Y1_nocap_l[r, :])
                Y2_l = p.capital_effect > 0.5 ? (@view p.basis_Y2_l[r, :]) : (@view p.basis_Y2_nocap_l[r, :])
                A_g, C_g = closed_form_ols(th_g, Y1_g, Y2_g, (@view p.w_g[r, :]), p.basis_C_max[r, 1])
                A_l, C_l = closed_form_ols(th_l, Y1_l, Y2_l, (@view p.w_l[r, :]), p.basis_C_max[r, 2])

                i = 0
                for h in d.horizon_dim
                    i += 1
                    zeta_g[i] = A_g * p.basis_Z1_g[r, h] + C_g * p.basis_Z2_g[r, h]
                    zeta_l[i] = A_l * p.basis_Z1_l[r, h] + C_l * p.basis_Z2_l[r, h]
                end

                z_g = 0.0
                z_l = 0.0
                for s in 0:min(10, t.t - start_idx)
                    idx = TimestepIndex(t.t - s)
                    # 全球渠道 = 受迫全球温度偏离 + 全球年际变率。回归端的
                    # g_shock 就是全球温度的 Hamilton 残差（即变率），而 §6.5
                    # 把局地温度中随其同向变动的部分也归给了这条渠道，故 gv
                    # 必须在这里、而不是在 z_l 里（open_issues.md P1-19）。
                    # 注意喂的是 gv 本身而非 coef_gvtas·gv：各国对全球变率的
                    # 放大已隐含在逐国的 zeta_g 里，与 T_global 的处理一致。
                    z_g += zeta_g[s + 1] * (p.T_global[idx] - T0 + p.global_lv[idx])
                    z_l += zeta_l[s + 1] * p.region_lv_local[idx, r]
                end
                v.z_tfp[t, r] = z_g + z_l

                if p.capital_effect > 0.5
                    rho = p.rho_c[r]
                    omega = 1.0 - rho
                    if t.t > start_idx
                        prev_t = TimestepIndex(t.t - 1)
                        v.khat[t, r] = rho * v.khat[prev_t, r] + omega * v.yhat[prev_t, r]
                    else
                        v.khat[t, r] = 0.0
                    end
                    v.yhat[t, r] = v.z_tfp[t, r] + p.alpha_c[r] * v.khat[t, r]
                else
                    v.khat[t, r] = 0.0
                    v.yhat[t, r] = v.z_tfp[t, r]
                end
                raw_ypc = p.ypc_baseline[t, r] * exp(v.yhat[t, r])
                safe_max = max(200.1, p.ypc_baseline[t, r] * 2.0)
                v.clipped[t, r] = (raw_ypc < 200.0 || raw_ypc > safe_max) ? 1.0 : 0.0
                v.ypc_new[t, r] = clamp(raw_ypc, 200.0, safe_max)
                v.loss[t, r] = (p.ypc_baseline[t, r] - v.ypc_new[t, r]) * p.population[t, r]
            end
        end
    end
end
