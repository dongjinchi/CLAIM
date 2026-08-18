# ============================================================================
# 【climateregional 组件】把全球平均温度转换为国家温度。
#
# 分段构造（2026-07-30 修改，open_issues.md P1-16）
# ---------------------------------------------------------------------------
# 历史段（year <= observed_end_year，生产设定为 2024）**直接使用观测**：
#
#     regtmp = regtmp_observed          （Berkeley 人口加权国家温度，LP 同源）
#
# 此前历史段也走 MESMER 重构（β·(T_FaIR − T_1950) + γ + 1950 单年锚点），
# 而观测数据本来就存在。实测该重构相对观测的误差：逐国逐年 sd 0.561 °C、
# |max| 3.18 °C，89/185 国的 RMSE 超过 0.5 °C；驱动损害起始年 T_state 的
# 1995–2024 气候态偏差 |max| 1.47 °C，52 国超过 0.5 °C。系统性偏冷 −0.18 °C
# 来自 FaIR 历史全球路径略低于 Berkeley 观测，以及 MESMER 的 γ 截距。
#
# 前向段（year > observed_end_year）用 MESMER 斜率 + **观测锚定的水平**：
#
#     regtmp = β_c · (T_global − T̄_global) + ō_c + region_lv
#
# 其中 T̄_global 与 ō_c 分别是 FaIR 全球温度和观测国家温度在锚定窗口
# （拼接年往前 anchor_window 年）上的均值。
#
# 为什么锚在 30 年均值而不是拼接年单年值
# ---------------------------------------------------------------------------
# 观测 = 趋势 + 天气。拼接年单年值含该年的天气异常，若用它做锚，这个异常会被
# 固化进前向趋势线的水平、永不衰减，而模型还会另外抽 ν 加上去——等于把同一年
# 的天气算两遍。实测 2024 是破纪录暖年，高出趋势 +0.368 °C，单年锚会让整条
# 2025–2150 路径永久偏高这么多。
#
# 30 年均值锚把天气平均掉。注意 water level 与参考点必须取**同一纪元**：
# ō_c 取观测窗口均值时，T̄_global 也必须取同窗口的 FaIR 均值，两者的纪元偏移
# 才会相互抵消。该式与"固定斜率 β 的 OLS 截距拟合"在代数上完全等价
# （实测最大差 1.8e-14）。
#
# 由此拼接处相对**观测值**会有约 −0.34 °C 的落差，这是刻意的——它等于减去
# 拼接年的天气异常。观测数据里 41.9% 的相邻年份变化都超过 0.34 °C，故该落差
# 在气候序列中稀松平常；而且前向段还会叠加 ν（国家级 sd 约 0.34），模拟的
# 首年并非确定地低于观测末年。要接续的是**趋势**，不是实现值。
#
# 连带后果：MESMER 的 region_lt_inter（γ）与 tempbase 在前向段不再被使用，
# 水平改由观测拟合——这是标准的 delta-change 偏差校正，保留 MESMER 的敏感度
# 模式（β）、用观测校正水平。两者仍保留为诊断输出与 legacy 分支所用。
#
# 无观测数据时（observed_end_year 早于模型起始年，合成数据测试即如此）退回
# legacy 构造：anchor 取 T_global_ref 与 tempbase + region_lt_inter。该分支由
# observed_end_year 显式选择，不是静默 fallback。
# ============================================================================
@defcomp climateregional begin
    regions   = Index()

    region_lt_coef  = Parameter(index = [regions])
    region_lt_inter = Parameter(index = [regions])
    region_lv   = Parameter(index = [time,regions])
    inputtemp = Parameter(index=[time])
    T_global_ref = Parameter()
    tempbase = Parameter(index=[regions])

    # 观测国家温度（Berkeley 人口加权，与 LP 的 lctmp_bkly_pw 同源）。
    # 仅 year <= observed_end_year 的行有效，其余年份的取值不被读取。
    regtmp_observed = Parameter(index = [time,regions])
    model_start_year = Parameter()
    observed_end_year = Parameter()
    anchor_window = Parameter()

    ts_lt = Variable(index = [time,regions])
    ts = Variable(index = [time,regions])
    regtmp = Variable(index=[time,regions])
    # 前向段的水平锚：观测窗口均值 ō_c。在拼接年算定后逐年沿用。
    anchor_base = Variable(index=[regions])
    anchor_tglobal = Variable(index=[regions])

    function run_timestep(p, v, d, t)
        year = Int(p.model_start_year) + t.t - 1
        splice = Int(p.observed_end_year)
        window = Int(p.anchor_window)

        # 锚定窗口 [splice-window+1, splice] 在模型时间轴上的下标。
        # 窗口左端可能早于模型起始年，下面统一用 max(1, lo) 截断为扩张窗口。
        lo = splice - window + 1 - Int(p.model_start_year) + 1
        hi = splice - Int(p.model_start_year) + 1
        # 只要拼接年不早于模型起始年，就存在可用的观测段。
        has_observed = hi >= 1

        if year <= splice && has_observed
            # ---- 历史段：直接用观测，不重构 ----
            for r in d.regions
                v.regtmp[t, r] = p.regtmp_observed[t, r]
                # ts_lt / ts 仅作诊断：观测已含天气，无法再拆分，故记为
                # 相对 tempbase 的偏离与其自身。
                v.ts_lt[t, r] = v.regtmp[t, r] - p.tempbase[r]
                v.ts[t, r] = v.ts_lt[t, r]
            end

            # 走到拼接年（或模型末年早于拼接年时的最后一年）就算定锚。
            # 每个历史时间步都重算一次，最终留下的即所需值。
            if t.t >= max(1, lo)
                tg = 0.0
                n = 0
                for i in max(1, lo):min(hi, t.t)
                    tg += p.inputtemp[TimestepIndex(i)]
                    n += 1
                end
                tg /= n
                for r in d.regions
                    acc = 0.0
                    for i in max(1, lo):min(hi, t.t)
                        acc += p.regtmp_observed[TimestepIndex(i), r]
                    end
                    v.anchor_base[r] = acc / n
                    v.anchor_tglobal[r] = tg
                end
            end
        else
            # ---- 前向段：MESMER 斜率 + 观测锚定的水平 ----
            for r in d.regions
                base, ref = if has_observed
                    v.anchor_base[r], v.anchor_tglobal[r]
                else
                    # legacy 构造（合成数据测试）：1950 单年锚 + MESMER 截距
                    p.tempbase[r] + p.region_lt_inter[r], p.T_global_ref
                end
                v.ts_lt[t, r] = p.region_lt_coef[r] * (p.inputtemp[t] - ref) + base - p.tempbase[r]
                v.ts[t, r] = v.ts_lt[t, r] + p.region_lv[t, r]
                v.regtmp[t, r] = v.ts[t, r] + p.tempbase[r]
            end
        end
    end
end
