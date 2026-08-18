# ============================================================================
# 【bhm_damage 组件】Burke/Hsiang/Miguel 二次温度-增长损害组件。
# 使用绝对温度而非 anomaly/shock：年度增长损害为
# f(T_climate_abs) - f(T_baseline_abs)，其中 baseline_abs = tempbase + region_lv。
# 当前版本不对 BHM 系数做不确定性抽样。
# ============================================================================
@defcomp bhm_damage begin
    regions = Index()

    regtmp = Parameter(index=[time, regions])
    region_lv = Parameter(index=[time, regions])
    tempbase = Parameter(index=[regions])
    ypc_baseline = Parameter(index=[time, regions])
    population = Parameter(index=[time, regions])
    beta1 = Parameter()
    beta2 = Parameter()
    model_start_year = Parameter()
    damage_start_year = Parameter()

    khat = Variable(index=[time, regions])
    ypc_new = Variable(index=[time, regions])
    loss = Variable(index=[time, regions])
    yhat = Variable(index=[time, regions])
    z_tfp = Variable(index=[time, regions])
    clipped = Variable(index=[time, regions])

    function run_timestep(p, v, d, t)
        start_idx = max(1, Int(p.damage_start_year) - Int(p.model_start_year) + 1)
        for r in d.regions
            if t.t < start_idx
                v.khat[t, r] = 0.0
                v.ypc_new[t, r] = p.ypc_baseline[t, r]
                v.loss[t, r] = 0.0
                v.yhat[t, r] = 0.0
                v.z_tfp[t, r] = 0.0
                v.clipped[t, r] = 0.0
            else
                T_climate = p.regtmp[t, r]
                T_baseline = p.tempbase[r] + p.region_lv[t, r]
                annual_damage = p.beta1 * (T_climate - T_baseline) + p.beta2 * (T_climate^2 - T_baseline^2)
                if t.t > start_idx
                    prev_t = TimestepIndex(t.t - 1)
                    v.z_tfp[t, r] = v.z_tfp[prev_t, r] + annual_damage
                else
                    v.z_tfp[t, r] = annual_damage
                end
                v.khat[t, r] = 0.0
                v.yhat[t, r] = v.z_tfp[t, r]
                raw_ypc = p.ypc_baseline[t, r] * exp(v.yhat[t, r])
                safe_max = max(200.1, p.ypc_baseline[t, r] * 2.0)
                v.clipped[t, r] = (raw_ypc < 200.0 || raw_ypc > safe_max) ? 1.0 : 0.0
                v.ypc_new[t, r] = clamp(raw_ypc, 200.0, safe_max)
                v.loss[t, r] = (p.ypc_baseline[t, r] - v.ypc_new[t, r]) * p.population[t, r]
            end
        end
    end
end
