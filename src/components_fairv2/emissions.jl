@defcomp emissions begin
    flourinated_gases  = Index()
    montreal_gases  = Index()
    aerosol_plus_gases  = Index()
    
    e_co2          = Parameter(index=[time])   # Annual carbon dioxide emissions (GtC yr⁻¹).
    e_ch4          = Parameter(index=[time])
    e_n2o          = Parameter(index=[time])
    e_flourinated  = Parameter(index=[time,flourinated_gases])
    e_montreal     = Parameter(index=[time,montreal_gases])
    e_aerosol_plus = Parameter(index=[time,aerosol_plus_gases])
    
    E_co2          = Variable(index=[time])
    E_ch4          = Variable(index=[time])
    E_n2o          = Variable(index=[time])
    E_flourinated  = Variable(index=[time,flourinated_gases])
    E_montreal     = Variable(index=[time,montreal_gases])
    E_aerosol_plus = Variable(index=[time,aerosol_plus_gases])
    
    function run_timestep(p, v, d, t)
        
        v.E_co2[t] = p.e_co2[t]
        v.E_ch4[t] = p.e_ch4[t]
        v.E_n2o[t] = p.e_n2o[t]
        
        for g in d.flourinated_gases
            v.E_flourinated[t,g] = p.e_flourinated[t,g]
        end
        
        for g in d.montreal_gases
            v.E_montreal[t,g] = p.e_montreal[t,g]
        end
        
        for g in d.aerosol_plus_gases
             v.E_aerosol_plus[t,g] = p.e_aerosol_plus[t,g]
        end
    end
end
