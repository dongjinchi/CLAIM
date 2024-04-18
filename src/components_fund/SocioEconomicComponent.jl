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
                #println(t,",",r,": ", v.ypc[t,r],v.ypc[t-1,r], v.ypcgrowth[t, r], "......")
            end
        end
    end
end