# ============================================================================
# 【socioeconomic 组件】外生社会经济基准模块。
# 由初始人均收入 ypc0 与外生增长率路径 ypcgrowth 递推每期人均收入：
#   ypc[t] = ypc[t-1] * (1 + ypcgrowth[t])，income = ypc * population。
# 基准路径刻意不受气候影响——气候损害由下游 claim_damage 组件以"偏离基准"的方式体现。
# ============================================================================
@defcomp socioeconomic begin
	regions = Index()

    income = Variable(index=[time,regions])
    ypc    = Variable(index=[time,regions]) 
    
    ypc0      = Parameter(index=[regions])
    ypcgrowth = Parameter(index=[time,regions])
	population = Parameter(index=[time,regions])

    function run_timestep(p, v, d, t)
        
        for r in d.regions
            
            if is_first(t)
                v.ypc[t,r] = p.ypc0[r]
                v.income[t,r] = p.ypc0[r] * p.population[t,r]
            else
                v.ypc[t,r] = v.ypc[t-1,r] * (1+p.ypcgrowth[t,r])
                v.income[t,r] = v.ypc[t,r] * p.population[t,r]
            end
        end
    end
end