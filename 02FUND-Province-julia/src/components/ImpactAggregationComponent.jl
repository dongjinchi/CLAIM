@defcomp impactaggregation begin
    regions = Index()
    
    scengrowth1 = Variable(index = [time, regions])
    #scengrowthchange = Variable(index = [time, regions])
    lnloss = Variable(index = [time, regions])
    loss = Variable(index = [time, regions])
    
    regtmp = Parameter(index=[time, regions])
    ypcgrowth = Parameter(index=[time, regions])
    income = Parameter(index=[time, regions])
    linearterm = Parameter()
    quadrterm = Parameter()
    ypcgrowth0 = Parameter(index=[regions])
    
    function run_timestep(p,v,d,t)
            
        for r in d.regions
            if is_first(t)
                p.ypcgrowth[t,r] = p.ypcgrowth0[r]
                v.scengrowth1[t,r] = p.linearterm + 2 * p.quadrterm * p.regtmp[t,r]
                #v.scengrowthchange[t,r] = v.scengrowth[t,r] - p.growth0[r]
                # Calculate the losses in logarithmic form
                v.lnloss[t,r] = log(1 + 0.01 * p.ypcgrowth[t,r] + v.scengrowth1[t,r])-log(1 + 0.01 * p.ypcgrowth[t,r])
                v.loss[t,r] = (1-exp(v.lnloss[t,r])) * p.income[t,r] * 10^9
            else
                v.scengrowth1[t,r] = p.linearterm + 2 * p.quadrterm * p.regtmp[t,r]
                #v.scengrowthchange[t,r] = v.scengrowth[t,r] - v.scengrowth[t-1,r]
                # Calculate the losses in logarithmic form
                v.lnloss[t,r] = v.lnloss[t-1,r] + log(1+ 0.01 * p.ypcgrowth[t,r] + v.scengrowth1[t,r]) - log(1 + 0.01 * p.ypcgrowth[t,r])
                v.loss[t,r] = (1 - exp(v.lnloss[t,r])) * p.income[t,r] * 10^9
            end
        end
    end
end
            
