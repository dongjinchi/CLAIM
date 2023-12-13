@defcomp climateregional begin
    regions = Index()
    inputtemp = Parameter(index=[time])
    bregtmp = Parameter(default = 1)
    #bregstmp = Parameter(index=[regions])
    temp0 = Parameter(index=[regions])

    temp = Variable(index=[time,regions]) #global temperautre mean tempearture above pre-industrial (it is same across regions)
    regtmp = Variable(index=[time,regions]) #regional temperautre
    regstmp = Variable(index=[time,regions]) 

    function run_timestep(p, v, d, t)

        for r in d.regions
            v.regtmp[t, r] = p.inputtemp[t] * p.bregtmp +p.temp0[r]
        end

        for r in d.regions
            v.temp[t, r] = v.regtmp[t, r] / p.bregtmp
        end

        #for r in d.regions
        #    v.regstmp[t, r] = p.inputtemp[t] * p.bregstmp[r] 
        #end
    end
end
