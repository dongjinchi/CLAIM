@defcomp impactaggregation begin
    regions = Index()
    
    growthchange = Variable(index = [time, regions])
    #scengrowthchange = Variable(index = [time, regions])
    lnloss = Variable(index = [time, regions])
    loss = Variable(index = [time, regions])
    ypc1 = Variable(index = [time, regions])
    temp2019 = Parameter(index=[regions])
    
    regtmp = Parameter(index=[time, regions])
    ypcgrowth = Parameter(index=[time, regions])
    income = Parameter(index=[time, regions])
    linearterm = Parameter()
    quadrterm = Parameter()
    ypc = Parameter(index=[time, regions])
    #ypcgrowth0 = Parameter(index=[regions])
    
    function run_timestep(p,v,d,t)
            
        for r in d.regions
            if t < TimestepValue(2020)
                v.lnloss[t,r] = 0
                v.growthchange[t,r] =0
                v.ypc1[t,r] = p.ypc[t,r]
            else
                v.growthchange[t,r] = p.linearterm*(p.regtmp[t,r] - p.temp2019[r]) + p.quadrterm*(p.regtmp[t,r]^2 - p.temp2019[r]^2)
                #Check for unrealistic values
                if v.ypc1[t-1,r] >= 10^6
                    v.growthchange[t,r] = 0
                end
                
                v.ypc1[t,r] = v.ypc1[t-1,r] * (1+ p.ypcgrowth[t,r] + v.growthchange[t,r])
                
                    
                v.lnloss[t,r] = v.lnloss[t-1,r] + log(1+ p.ypcgrowth[t,r] + v.growthchange[t,r]) - log(1 + p.ypcgrowth[t,r])
                v.loss[t,r] = (1 - exp(v.lnloss[t,r])) * p.income[t,r]
            end
        end
    end
end