@defcomp climateregional begin
    regions   = Index()
    
    region_lt_coef  = Parameter(index = [regions])
    region_lt_inter = Parameter(index = [regions])
    region_lv   = Parameter(index = [time,regions])
    inputtemp = Parameter(index=[time])
    temp1990     = Parameter(index=[regions])
    
    ts_lt = Variable(index = [time,regions])
    ts = Variable(index = [time,regions])
    regtmp = Variable(index=[time,regions])
    temp2019 = Variable(index=[regions])

    function run_timestep(p, v, d, t)
        
        for r in d.regions
            if is_first(t)
                v.regtmp[t, r] = p.temp1990[r]
            else
                v.ts_lt[t,r] = p.region_lt_coef[r] * (p.inputtemp[t]-0.594499) + p.region_lt_inter[r]
                v.ts[t,r] = v.ts_lt[t,r]+p.region_lv[t,r]
                v.regtmp[t, r] = v.ts[t,r] +p.temp1990[r]
            end
            
            if t == TimestepValue(2019)
                v.temp2019[r] = v.regtmp[t,r]
            end
        end
    end
end
        
        
