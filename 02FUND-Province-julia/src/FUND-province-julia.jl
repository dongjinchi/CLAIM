module FUNDprovince

using Mimi
using DelimitedFiles
using Random

include("helper.jl")

include("new_marginaldamages.jl")
include("montecarlo/defmcs.jl")
include("montecarlo/run_fund_mcs.jl")
include("montecarlo/run_fund_scc_mcs.jl")

include("components/SocioEconomicComponent.jl")
include("components/PopulationComponent.jl")
include("components/EmissionsComponent.jl")
#include("components/GeographyComponent.jl")
include("components/ScenarioUncertaintyComponent.jl")
include("components/ClimateCO2CycleComponent.jl")
#include("components/ClimateCH4CycleComponent.jl")
#include("components/ClimateN2OCycleComponent.jl")
#include("components/ClimateSF6CycleComponent.jl")
include("components/ClimateForcingComponent.jl")
include("components/ClimateDynamicsComponent.jl")
#include("components/BioDiversityComponent.jl")
include("components/ClimateRegionalComponent.jl")
#include("components/OceanComponent.jl")
#include("components/ImpactCardiovascularRespiratoryComponent.jl")
#include("components/ImpactDiarrhoeaComponent.jl")
#include("components/ImpactExtratropicalStormsComponent.jl")
#include("components/ImpactDeathMorbidityComponent.jl")
#include("components/ImpactVectorBorneDiseasesComponent.jl")
#include("components/ImpactTropicalStormsComponent.jl")
#include("components/ImpactSeaLevelRiseComponent.jl")
include("components/ImpactAggregationComponent.jl") #replaced with a new aggregate impact component
#include("components/VslVmorbComponent.jl")


function get_model(;syear,  nsteps::Int, regions, datadir, params)
    
    # ------------------------------------
    # Create model
    # ------------------------------------
    
    m = Model()
    
    # ------------------------------------
    # Set dimensions
    # ------------------------------------
    
    set_dimension!(m, :time, collect(syear:syear+nsteps)) # set time dimension
    set_dimension!(m, :regions, regions)
    
    # ------------------------------------
    # Create components
    # ------------------------------------
    
    add_comp!(m, scenariouncertainty)
    add_comp!(m, population)
    #add_comp!(m, geography)
    add_comp!(m, socioeconomic)
    add_comp!(m, emissions)
    add_comp!(m, climateco2cycle)
    #add_comp!(m, climatech4cycle)
    #add_comp!(m, climaten2ocycle)
    #add_comp!(m, climatesf6cycle)
    add_comp!(m, climateforcing)
    add_comp!(m, climatedynamics)
    #add_comp!(m, biodiversity)
    add_comp!(m, climateregional)
    #add_comp!(m, ocean)
    #add_comp!(m, impactcardiovascularrespiratory)
    #add_comp!(m, impactdiarrhoea)
    #add_comp!(m, impactextratropicalstorms)
    #add_comp!(m, impactvectorbornediseases)
    #add_comp!(m, impacttropicalstorms)
    #add_comp!(m, vslvmorb)
    #add_comp!(m, impactdeathmorbidity)
    #add_comp!(m, impactsealevelrise)
    add_comp!(m, impactaggregation)

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    #connect_param!(m, :geography, :landloss, :impactsealevelrise, :landloss)

    connect_param!(m, :population, :pgrowth, :scenariouncertainty, :pgrowth)
    #connect_param!(m, :population, :enter, :impactsealevelrise, :enter)
    #connect_param!(m, :population, :leave, :impactsealevelrise, :leave)
    #connect_param!(m, :population, :dead, :impactdeathmorbidity, :dead)

    #connect_param!(m, :socioeconomic, :area, :geography, :area)
    connect_param!(m, :socioeconomic, :globalpopulation, :population, :globalpopulation)
    connect_param!(m, :socioeconomic, :populationin1, :population, :populationin1)
    connect_param!(m, :socioeconomic, :population, :population, :population)
    connect_param!(m, :socioeconomic, :pgrowth, :scenariouncertainty, :pgrowth)
    #connect_param!(m, :socioeconomic, :ypcgrowth, :scenariouncertainty, :ypcgrowth)
    connect_param!(m, :socioeconomic, :ygrowth, :scenariouncertainty, :ygrowth)
    connect_param!(m, :socioeconomic, :mitigationcost, :emissions, :mitigationcost)

    connect_param!(m, :emissions, :income, :socioeconomic, :income)
    connect_param!(m, :emissions, :population, :population, :population)
    #connect_param!(m, :emissions, :forestemm, :scenariouncertainty, :forestemm)
    connect_param!(m, :emissions, :aeei, :scenariouncertainty, :aeei)
    connect_param!(m, :emissions, :acei, :scenariouncertainty, :acei)
    #connect_param!(m, :emissions, :ypcgrowth, :scenariouncertainty, :ypcgrowth)

    connect_param!(m, :climateco2cycle, :mco2, :emissions, :mco2)

    #connect_param!(m, :climatech4cycle, :globch4, :emissions, :globch4)

    #connect_param!(m, :climaten2ocycle, :globn2o, :emissions, :globn2o)

    connect_param!(m, :climateco2cycle, :temp, :climatedynamics, :temp)

    #connect_param!(m, :climatesf6cycle, :globsf6, :emissions, :globsf6)

    connect_param!(m, :climateforcing, :acco2, :climateco2cycle, :acco2)
    #connect_param!(m, :climateforcing, :acch4, :climatech4cycle, :acch4)
    #connect_param!(m, :climateforcing, :acn2o, :climaten2ocycle, :acn2o)
    #connect_param!(m, :climateforcing, :acsf6, :climatesf6cycle, :acsf6)

    connect_param!(m, :climatedynamics, :radforc, :climateforcing, :radforc)

    connect_param!(m, :climateregional, :inputtemp, :climatedynamics, :temp)

    #connect_param!(m, :biodiversity, :temp, :climatedynamics, :temp)

    #connect_param!(m, :ocean, :temp, :climatedynamics, :temp)

    #connect_param!(m, :impactcardiovascularrespiratory, :population, :population, :population)
    #connect_param!(m, :impactcardiovascularrespiratory, :temp, :climateregional, :temp)
    #connect_param!(m, :impactcardiovascularrespiratory, :plus, :socioeconomic, :plus)
    #connect_param!(m, :impactcardiovascularrespiratory, :urbpop, :socioeconomic, :urbpop)

    #connect_param!(m, :impactdiarrhoea, :population, :population, :population)
    #connect_param!(m, :impactdiarrhoea, :income, :socioeconomic, :income)
    #connect_param!(m, :impactdiarrhoea, :regtmp, :climateregional, :regtmp)

    #connect_param!(m, :impactextratropicalstorms, :population, :population, :population)
    #connect_param!(m, :impactextratropicalstorms, :income, :socioeconomic, :income)
    #connect_param!(m, :impactextratropicalstorms, :acco2, :climateco2cycle, :acco2)

    #connect_param!(m, :impactvectorbornediseases, :population, :population, :population)
    #connect_param!(m, :impactvectorbornediseases, :income, :socioeconomic, :income)
    #connect_param!(m, :impactvectorbornediseases, :temp, :climateregional, :temp)

    #connect_param!(m, :impacttropicalstorms, :population, :population, :population)
    #connect_param!(m, :impacttropicalstorms, :income, :socioeconomic, :income)
    #connect_param!(m, :impacttropicalstorms, :regstmp, :climateregional, :regstmp)

    #connect_param!(m, :vslvmorb, :population, :population, :population)
    #connect_param!(m, :vslvmorb, :income, :socioeconomic, :income)

    #connect_param!(m, :impactdeathmorbidity, :vsl, :vslvmorb, :vsl)
    #connect_param!(m, :impactdeathmorbidity, :vmorb, :vslvmorb, :vmorb)
    #connect_param!(m, :impactdeathmorbidity, :population, :population, :population)
    #connect_param!(m, :impactdeathmorbidity, :dengue, :impactvectorbornediseases, :dengue)
    #connect_param!(m, :impactdeathmorbidity, :schisto, :impactvectorbornediseases, :schisto)
    #connect_param!(m, :impactdeathmorbidity, :malaria, :impactvectorbornediseases, :malaria)
    #connect_param!(m, :impactdeathmorbidity, :cardheat, :impactcardiovascularrespiratory, :cardheat)
    #connect_param!(m, :impactdeathmorbidity, :cardcold, :impactcardiovascularrespiratory, :cardcold)
    #connect_param!(m, :impactdeathmorbidity, :resp, :impactcardiovascularrespiratory, :resp)
    #connect_param!(m, :impactdeathmorbidity, :diadead, :impactdiarrhoea, :diadead)
    #connect_param!(m, :impactdeathmorbidity, :diasick, :impactdiarrhoea, :diasick)
    #connect_param!(m, :impactdeathmorbidity, :hurrdead, :impacttropicalstorms, :hurrdead)
    #connect_param!(m, :impactdeathmorbidity, :extratropicalstormsdead, :impactextratropicalstorms, :extratropicalstormsdead)

    #connect_param!(m, :impactsealevelrise, :population, :population, :population)
    #connect_param!(m, :impactsealevelrise, :income, :socioeconomic, :income)
    #connect_param!(m, :impactsealevelrise, :sea, :ocean, :sea)
    #connect_param!(m, :impactsealevelrise, :area, :geography, :area)
   
    connect_param!(m, :impactaggregation, :regtmp,:climateregional, :regtmp)
    connect_param!(m, :impactaggregation, :ypcgrowth, :socioeconomic, :ypcgrowth)
    connect_param!(m, :impactaggregation, :income,:socioeconomic, :income)
    
    # ---------------------------------------------
    # Set all external parameter values
    # ---------------------------------------------

    parameters = params === nothing ? load_default_parameters(nsteps, datadir) : params
    
    # Set unshared parameters - name is a Tuple{Symbol, Symbol} of (component_name, param_name)
    for (name, value) in parameters[:unshared]
        show(name[1])
        update_param!(m, name[1], name[2], value)
    end

    update_param!(m, :impactaggregation,:linearterm,  0.0162)
    update_param!(m, :impactaggregation,:quadrterm,  -0.000903)

    # Set shared parameters - name is a Symbol representing the param_name, here
    # we will create a shared model parameter with the same name as the component
    # parameter and then connect our component parameters to this shared model parameter
    
    # * for convenience later, name shared model parameter same as the component 
    # parameters, but this is not required could give a unique name *
    #add_shared_param!(m, :ch4pre, parameters[:shared][:ch4pre])
    #connect_param!(m, :climateforcing, :ch4pre, :ch4pre)
    #connect_param!(m, :climatech4cycle, :ch4pre, :ch4pre)

    #add_shared_param!(m, :n2opre, parameters[:shared][:n2opre])
    #connect_param!(m, :climateforcing, :n2opre, :n2opre)
    #connect_param!(m, :climaten2ocycle, :n2opre, :n2opre)

    #add_shared_param!(m, :sf6pre, parameters[:shared][:sf6pre])
    #connect_param!(m, :climateforcing, :sf6pre, :sf6pre)
    #connect_param!(m, :climatesf6cycle, :sf6pre, :sf6pre)

    add_shared_param!(m, :co2pre, parameters[:shared][:co2pre])
    connect_param!(m, :climateforcing, :co2pre, :co2pre)
    #connect_param!(m, :impactextratropicalstorms, :co2pre, :co2pre)

    #add_shared_param!(m, :nospecbase, parameters[:shared][:nospecbase])
    #connect_param!(m, :biodiversity, :nospecbase, :nospecbase)

    #add_shared_param!(m, :dbsta, parameters[:shared][:dbsta])
    #connect_param!(m, :biodiversity, :dbsta, :dbsta)

    #add_shared_param!(m, :bregtmp, parameters[:shared][:bregtmp], dims=[:regions])
    #connect_param!(m, :climateregional, :bregtmp, :bregtmp)
    #connect_param!(m, :impactdiarrhoea, :bregtmp, :bregtmp)

    #add_shared_param!(m, :plus90, parameters[:shared][:plus90], dims=[:regions])
    #connect_param!(m, :impactcardiovascularrespiratory, :plus90, :plus90)
    #connect_param!(m, :socioeconomic, :plus90, :plus90)

    #add_shared_param!(m, :gdp90, parameters[:shared][:gdp90], dims=[:regions])
    #connect_param!(m, :emissions, :gdp90, :gdp90)
    #connect_param!(m, :impactdiarrhoea, :gdp90, :gdp90)
    #connect_param!(m, :impactextratropicalstorms, :gdp90, :gdp90)
    #connect_param!(m, :impacttropicalstorms, :gdp90, :gdp90)
    #connect_param!(m, :impactvectorbornediseases, :gdp90, :gdp90)
    #connect_param!(m, :socioeconomic, :gdp90, :gdp90)

    #add_shared_param!(m, :pop90, parameters[:shared][:pop90], dims=[:regions])
    #connect_param!(m, :emissions, :pop90, :pop90)
    #connect_param!(m, :impactdiarrhoea, :pop90, :pop90)
    #connect_param!(m, :impactextratropicalstorms, :pop90, :pop90)
    #connect_param!(m, :impacttropicalstorms, :pop90, :pop90)
    #connect_param!(m, :impactvectorbornediseases, :pop90, :pop90)
    #connect_param!(m, :socioeconomic, :pop90, :pop90)
    
    return m

end
end



