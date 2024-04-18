# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------
# This file contains functions and other snippets of code that are used in various calculations for MimiFAIRv2.
# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------


#######################################################################################################################
# CALCULATE CLIMATE & THERMAL RESPONSE PARAMETERS
#######################################################################################################################
# Description: From original Python FAIR code, "This function returns the FaIRv2.0.0-alpha default climate parameters.
#              Use the kwargs to specify pre-determined climate sensitivities. In both cases, the response timescales d1-3
#              (and the shortest-timescale coefficient, q1) are set to the central estimate of a CMIP6 inferred distribution
#              constrained with observational warming. The constraint does not significantly affect the central estimates of
#              the prior (ie. raw CMIP6 inference) distribution."
#
# Function Arguments:
#
#       TCR: Transient climate response (K).
#       RWF: Realized warming fraction (ratio of TCR/ECS).
#       F2x: Radiative forcing from a doubling of carbon dioxide concentrations (Wm⁻²).
#----------------------------------------------------------------------------------------------------------------------

"""
    get_thermal_parameter_defaults(;TCR::Float64=1.79, RWF::Float64=0.552, F2x::Float64=3.759)

Return a dataframe of the FAIRv2.0.0-alpha default climate parameters using the
optional keyword arguments for transient climate response (K) (TCR; defaults to 1.79), 
realized warming fraction in ratio of TCR/ECS (RWF; defaults to 0.552), and radiative 
forcing from a double of carbon dioxide concentrations (Wm⁻²) (F2x; defaults to 
3.759).

The response timescales d1-3 (and the shortest-timescale coefficient, q1) are set 
to the central estimate of a CMIP6 inferred distribution constrained with observational 
warming. The constraint does not significantly affect the central estimates of
the prior (ie. raw CMIP6 inference) distribution.
"""
function get_thermal_parameter_defaults(;TCR::Float64=1.79, RWF::Float64=0.552, F2x::Float64=3.759)

    # Set default values for d1, d2, d3, and q1 parameters.
    d1 = 0.903
    d2 = 7.92
    d3 = 355
    q1 = 0.180

    # Calculate equilibrium climate sensitivity based on user-specified TCR and RWF values.
    ECS = TCR/RWF

    # Intermediate calculations.
    v1 = (1-(d1/69.66) * (1-exp(-69.66/d1)))
    v2 = (1-(d2/69.66) * (1-exp(-69.66/d2)))
    v3 = (1-(d3/69.66) * (1-exp(-69.66/d3)))

    # Calculate equilibrium response of second (q2) and third (q3) thermal boxes.
    q3 = (((TCR/F2x) - q1*(v1-v2) - (ECS/F2x)*v2) / (v3-v2))
    q2 = (ECS/F2x - q1 -  q3)

    # Return thermal response parameters as a dataframe.
    df = DataFrame(d1=d1, d2=d2, d3=d3, q1=q1, q2=q2, q3=q3, v1=v1, v2=v2, v3=v3, tcr=TCR, rwf=RWF, F2x=F2x, ECS=ECS)

    return df
end



#######################################################################################################################
# CALCULATE CONSTANTS TO APPROXIMATE STATE-DEPENDENT TIMESCALE ADJUSTMENT FACTOR
#######################################################################################################################
# Description: This function calculates constants for estimating the state-dependent adjustment coefficient of a
#              reservior's lifetime (α) such that it approximates the Millar et al. (2017) numerical solution for a
#              iIRF100 carbon cycle parameterization at α=1.
#
# Function Arguments:
#
#       a: Fraction of emissions entering iᵗʰ atmospheric pool (by default, FAIR has four pools).
#       τ: Atmospheric lifetime of gas in iᵗʰ pool.
#----------------------------------------------------------------------------------------------------------------------

# Function version for a single gas (where a and τ are vectors).
"""
    calculate_g0_g1(a::Array{Float64,1}, τ::Array{Float64,1})

For a SINGLE gas, calculate and return constants for estimating the state-dependent 
adjustment coefficient of a reservior's lifetime (α) such that it approximates 
the Millar et al. (2017) numerical solution for a iIRF100 carbon cycle parameterization 
at α=1. The two required arguments are a, the fraction of emissions entering iᵗʰ 
atmospheric pool, and τ, the atmospheric lifetime of gas in iᵗʰ pool.
"""
function calculate_g0_g1(a::Array{Float64,1}, τ::Array{Float64,1})
    g1 = sum(a .* τ .* (1.0 .- (1.0 .+ 100.0 ./ τ) .* exp.(-100.0 ./ τ)))
    g0 = exp(-1.0 * sum(a .* τ .* (1.0 .- exp.(-100.0 ./ τ))) / g1)
    return g0, g1
end

# Function version for multiple gases (where a and τ are 2-d arrays).
"""
    calculate_g0_g1(a::Array{Float64,1}, τ::Array{Float64,1})

For MULITPLE gases, calculate and return constants for estimating the state-dependent 
adjustment coefficient of a reservior's lifetime (α) such that it approximates 
the Millar et al. (2017) numerical solution for a iIRF100 carbon cycle parameterization 
at α=1. The two required arguments are a, the fraction of emissions entering iᵗʰ 
atmospheric pool, and τ, the atmospheric lifetime of gas in iᵗʰ pool.
"""
function calculate_g0_g1(a::Array{Float64,2}, τ::Array{Float64,2})
    g1 = vec(sum(a .* τ .* (1.0 .- (1.0 .+ 100.0 ./ τ) .* exp.(-100.0 ./ τ)), dims=2))
    g0 = vec(exp.(-1.0 .* sum(a .* τ .* (1.0 .- exp.(-100.0 ./ τ)), dims=2) ./ g1))
    return g0, g1
end



#######################################################################################################################
# CALCULATE RADIATIVE FORCING
#######################################################################################################################
# Description: This function calculates the effective radiative forcing based on the atmospheric concentration of a 
#              greenhouse gas or other forcing agent. It follows Equation (6) in Leach et al. (2021).
#
# Function Arguments:
#
#       f:    A vector of three concentration-forcing coefficients (1 = logarithmic, 2 = linear, 3 = square root).
#       C:    Concentration of forcing agent in atmosphere.
#       C_pi: Pre-industrial concentration of forcing agent in atmosphere.
#----------------------------------------------------------------------------------------------------------------------

"""
    calculate_RF(f, C, C_pi)

Calculate and return the effective radiative forcing based on the atmospheric 
concentration of a greenhouse gas or other forcing agent following Equation (6) 
in Leach et al. (2021). 

The required arguments are as follows:
    f:    A vector of three concentration-forcing coefficients (1 = logarithmic, 2 = linear, 3 = square root)
    C:    Concentration of forcing agent in atmosphere
    C_pi: Pre-industrial concentration of forcing agent in atmosphere
"""
function calculate_RF(f, C, C_pi)

    if C <= 0.0
        # If concentration is negative, set log and square-root concentrations to 0.0. This follows the original Python FAIRv2.0 code.
        forcing = f[1] * 0.0 + f[2] * (C - C_pi) + f[3] * (0.0 - sqrt(C_pi))
    else
        # Otherwise, calculate forcing as the sum of logarithmic, linear, and square-root terms.
        forcing = f[1] * log(C / C_pi) + f[2] * (C - C_pi) + f[3] * (sqrt(C) - sqrt(C_pi))
    end

    return forcing
end

#######################################################################################################################
# UPDATING PARAMETERS FOR MIMIFAIR2
#######################################################################################################################
# Description: This function updates the parameters for MimiFAIR2
#
# Function Arguments:
#
#       f:   
#----------------------------------------------------------------------------------------------------------------------

"""
    update_MimiFAIR2_params!(m; ar6_scenario::String="ssp245", start_year::Int=1750, end_year::Int=2300)

Update all parameter settings in the model `m` with default FAIR  v2 settings 
as pulled from the MimiFAIRv2 repository.  Optional arguments include `ar6_scenario` to indicate
what scenario to use, `start_year` and `end_year` to indicate the model time dimension.
"""

using CSVFiles, DataFrames, Mimi

function update_MimiFAIR2_params!(m; emissions_forcing_scenario::String="ssp119", start_year::Int=1950, end_year::Int=2300, TCR::Float64=1.79, RWF::Float64=0.552, F2x::Float64=3.759)

 	# ---------------------------------------------
	# ---------------------------------------------
	# Set Up Data and Parameter Values
	# ---------------------------------------------
 	# ---------------------------------------------

	# Load emissions and forcing data and crop to appropriate model years (scenarios span 1750-2500 by default).
	scenario_indices = indexin(start_year:end_year, 1750:2500)
	forcing_data     = DataFrame(load(joinpath(@__DIR__, "..", "data/fair_data", "rcmip_"*emissions_forcing_scenario*"_effective_radiative_forcing_1750_to_2500.csv"), skiplines_begin=6))[scenario_indices,:]
	emissions_data   = DataFrame(load(joinpath(@__DIR__, "..", "data/fair_data", "rcmip_"*emissions_forcing_scenario*"_emissions_1750_to_2500_rebased.csv"), skiplines_begin=6))[scenario_indices,:]

	# Load FAIR default gas cycle (gas) and indirect radiative forcing (irf_p) parameters.
	gas_p = DataFrame(load(joinpath(@__DIR__, "..", "data/fair_data", "default_gas_cycle_parameters.csv"), skiplines_begin=6))
	irf_p = DataFrame(load(joinpath(@__DIR__, "..", "data/fair_data", "default_indirect_radiative_forcing_parameters.csv"), skiplines_begin=7))

	# Load initial conditions for 1950 under the RCMIP emissions & forcing scenario.
	init_gas_vals     = DataFrame(load(joinpath(@__DIR__, "..", "data/fair_data", "fair_initial_gas_cycle_conditions_1950.csv"), skiplines_begin=7))
	init_thermal_vals = DataFrame(load(joinpath(@__DIR__, "..", "data/fair_data", "fair_initial_thermal_conditions_1950.csv"), skiplines_begin=7))

	# Isolate specific gas parameters for convenience.
	co2_p              = filter(:gas_name  => ==("carbon_dioxide"), gas_p)
	ch4_p              = filter(:gas_name  => ==("methane"), gas_p)
	n2o_p              = filter(:gas_name  => ==("nitrous_oxide"), gas_p)
	montreal_gas_p     = filter(:gas_group => ==("montreal"), gas_p)
	flourinated_gas_p  = filter(:gas_group => ==("flourinated"), gas_p)
	aerosol_plus_gas_p = filter(:gas_group => ==("aerosol_plus"), gas_p)

	# Sort arrays of parameters for multiple gases so they are listed alphabetically.
	sort!(flourinated_gas_p,  :gas_name)
	sort!(montreal_gas_p,     :gas_name)
	sort!(aerosol_plus_gas_p, :gas_name)

	#Isolate Montreal indirect forcing parameters and sort alphabetically.
	montreal_irf_p = filter(:gas_group => ==("montreal"), irf_p)
	sort!(montreal_irf_p, :gas_name)

	# Isolate initial model conditions (initial conditions currently default to year 1750).
	co2_init          = filter(:gas_name  => ==("carbon_dioxide"), init_gas_vals)
	ch4_init          = filter(:gas_name  => ==("methane"), init_gas_vals)
	n2o_init          = filter(:gas_name  => ==("nitrous_oxide"), init_gas_vals)
	montreal_init     = filter(:gas_group => ==("montreal"), init_gas_vals)
	flourinated_init  = filter(:gas_group => ==("flourinated"), init_gas_vals)
	aerosol_plus_init = filter(:gas_group => ==("aerosol_plus"), init_gas_vals)

	# Sort arrays of initial conditions for multiple gases so they are listed alphabetically.
	sort!(montreal_init,     :gas_name)
	sort!(flourinated_init,  :gas_name)
	sort!(aerosol_plus_init, :gas_name)

	# Extract emissions arrays for multi-gas groupings.
	montreal_emissions     = emissions_data[:, Symbol.(montreal_init.gas_name)]
	flourinated_emissions  = emissions_data[:, Symbol.(flourinated_init.gas_name)]
	aerosol_plus_emissions = emissions_data[:, Symbol.(aerosol_plus_init.gas_name)]

    # Create helper arrays of indices to use for pulling parameter groups.
    a_idxs = [:a1,:a2,:a3,:a4]
    τ_idxs = [:tau1,:tau2,:tau3,:tau4]

	# Create arrays for 'a' parameter groups. 
	co2_a          = vec(Array(co2_p[:, a_idxs]))
	ch4_a          = vec(Array(ch4_p[:, a_idxs]))
	n2o_a          = vec(Array(n2o_p[:, a_idxs]))
	montreal_a     = Array(montreal_gas_p[:, a_idxs])
	flourinated_a  = Array(flourinated_gas_p[:, a_idxs])
	aerosol_plus_a = Array(aerosol_plus_gas_p[:, a_idxs])

	# Create arrays for 'τ' parameter groups.
	co2_τ 		   = vec(Array(co2_p[:, τ_idxs]))
	ch4_τ 		   = vec(Array(ch4_p[:, τ_idxs]))
	n2o_τ 		   = vec(Array(n2o_p[:, τ_idxs]))
	montreal_τ 	   = Array(montreal_gas_p[:, τ_idxs])
	flourinated_τ  = Array(flourinated_gas_p[:, τ_idxs])
	aerosol_plus_τ = Array(aerosol_plus_gas_p[:, τ_idxs])

	# Calculate constants to approximate numerical solution for state-dependent timescale adjustment factor from Millar et al. (2017).
	g0_co2, g1_co2 					 =  calculate_g0_g1(co2_a, co2_τ)
	g0_ch4, g1_ch4 				     =  calculate_g0_g1(ch4_a, ch4_τ)
	g0_n2o, g1_n2o 					 =  calculate_g0_g1(n2o_a, n2o_τ)
	g0_montreal, g1_montreal 		 =  calculate_g0_g1(montreal_a, montreal_τ)
	g0_flourinated, g1_flourinated   =  calculate_g0_g1(flourinated_a, flourinated_τ)
	g0_aerosol_plus, g1_aerosol_plus =  calculate_g0_g1(aerosol_plus_a, aerosol_plus_τ)

	# Calculate default thermal parameter values 
    # Defaults pulled from get_model function defaults are are as follows:
    #   transient climate response = 1.79
    #   realized warming fraction = 0.552
    #   forcing from a doubling of CO₂ = 3.759
	thermal_p = get_thermal_parameter_defaults(TCR = TCR, RWF= RWF, F2x = F2x)

	# Calculate thermal decay factors, defined as exp(-1/d).
	thermal_decay_factors = exp.(-1.0 ./ vec(Array(thermal_p[:,[:d1,:d2,:d3]])))

	# ---------------------------------------------
    # Set component-specific parameters
    # ---------------------------------------------
    
    # --- Emissions --- #
	update_param!(m, :emissions, :e_co2, emissions_data.carbon_dioxide) 
	update_param!(m, :emissions, :e_ch4, emissions_data.methane)
	update_param!(m, :emissions, :e_n2o, emissions_data.nitrous_oxide)
	update_param!(m, :emissions, :e_flourinated, Array(flourinated_emissions))
	update_param!(m, :emissions, :e_montreal, Array(montreal_emissions))
	update_param!(m, :emissions, :e_aerosol_plus, Array(aerosol_plus_emissions))

 	# ---- Carbon Cycle ---- #
	update_param!(m, :co2_cycle, :co2_0, co2_init.concentration[1])
	update_param!(m, :co2_cycle, :r0_co2, co2_p.r0[1])
	update_param!(m, :co2_cycle, :rU_co2, co2_p.rC[1])
	update_param!(m, :co2_cycle, :rT_co2, co2_p.rT[1])
	update_param!(m, :co2_cycle, :rA_co2, co2_p.rA[1])
	update_param!(m, :co2_cycle, :GU_co2_0, co2_init.cumulative_uptake[1])
	update_param!(m, :co2_cycle, :g0_co2, g0_co2)
	update_param!(m, :co2_cycle, :g1_co2, g1_co2)
	update_param!(m, :co2_cycle, :R0_co2, vec(Array(co2_init[:,[:pool1,:pool2,:pool3,:pool4]])))
	update_param!(m, :co2_cycle, :emiss2conc_co2, co2_p.emis2conc[1])
	update_param!(m, :co2_cycle, :a_co2, co2_a)
	update_param!(m, :co2_cycle, :τ_co2, co2_τ)


	# ---- Methane Cycle ---- #
	update_param!(m, :ch4_cycle, :ch4_0, ch4_init.concentration[1])
	update_param!(m, :ch4_cycle, :r0_ch4, ch4_p.r0[1])
	update_param!(m, :ch4_cycle, :rU_ch4, ch4_p.rC[1])
	update_param!(m, :ch4_cycle, :rT_ch4, ch4_p.rT[1])
	update_param!(m, :ch4_cycle, :rA_ch4, ch4_p.rA[1])
	update_param!(m, :ch4_cycle, :GU_ch4_0, ch4_init.cumulative_uptake[1])
	update_param!(m, :ch4_cycle, :g0_ch4, g0_ch4)
	update_param!(m, :ch4_cycle, :g1_ch4, g1_ch4)
	update_param!(m, :ch4_cycle, :R0_ch4, vec(Array(ch4_init[:,[:pool1,:pool2,:pool3,:pool4]])))
	update_param!(m, :ch4_cycle, :emiss2conc_ch4, ch4_p.emis2conc[1])
	update_param!(m, :ch4_cycle, :a_ch4, ch4_a)
	update_param!(m, :ch4_cycle, :τ_ch4, ch4_τ)

	# ---- Nitrous Oxide Cycle ---- #
	update_param!(m, :n2o_cycle, :n2o_0, n2o_init.concentration[1])
	update_param!(m, :n2o_cycle, :r0_n2o, n2o_p.r0[1])
	update_param!(m, :n2o_cycle, :rU_n2o, n2o_p.rC[1])
	update_param!(m, :n2o_cycle, :rT_n2o, n2o_p.rT[1])
	update_param!(m, :n2o_cycle, :rA_n2o, n2o_p.rA[1])
	update_param!(m, :n2o_cycle, :GU_n2o_0, n2o_init.cumulative_uptake[1])
	update_param!(m, :n2o_cycle, :g0_n2o, g0_n2o)
	update_param!(m, :n2o_cycle, :g1_n2o, g1_n2o)
	update_param!(m, :n2o_cycle, :R0_n2o, vec(Array(n2o_init[:,[:pool1,:pool2,:pool3,:pool4]])))
	update_param!(m, :n2o_cycle, :emiss2conc_n2o, n2o_p.emis2conc[1])
	update_param!(m, :n2o_cycle, :a_n2o, n2o_a)
	update_param!(m, :n2o_cycle, :τ_n2o, n2o_τ)

	# ---- Flourinated Gas Cycles ---- #
	update_param!(m, :flourinated_cycles, :flourinated_0, flourinated_init[:, :concentration])
	update_param!(m, :flourinated_cycles, :r0_flourinated, flourinated_gas_p[:,:r0])
	update_param!(m, :flourinated_cycles, :rU_flourinated, flourinated_gas_p[:,:rC])
	update_param!(m, :flourinated_cycles, :rT_flourinated, flourinated_gas_p[:,:rT])
	update_param!(m, :flourinated_cycles, :rA_flourinated, flourinated_gas_p[:,:rA])
	update_param!(m, :flourinated_cycles, :GU_flourinated_0, flourinated_init[:, :cumulative_uptake])
	update_param!(m, :flourinated_cycles, :g0_flourinated, g0_flourinated)
	update_param!(m, :flourinated_cycles, :g1_flourinated, g1_flourinated)
	update_param!(m, :flourinated_cycles, :R0_flourinated, Array(flourinated_init[:,[:pool1,:pool2,:pool3,:pool4]]))
	update_param!(m, :flourinated_cycles, :emiss2conc_flourinated, flourinated_gas_p[:,:emis2conc])
	update_param!(m, :flourinated_cycles, :a_flourinated, flourinated_a)
	update_param!(m, :flourinated_cycles, :τ_flourinated, flourinated_τ)

	# ---- Montreal Protocol Gas Cycles ---- #
	update_param!(m, :montreal_cycles, :montreal_0, montreal_init[:, :concentration])
	update_param!(m, :montreal_cycles, :r0_montreal, montreal_gas_p[:,:r0])
	update_param!(m, :montreal_cycles, :rU_montreal, montreal_gas_p[:,:rC])
	update_param!(m, :montreal_cycles, :rT_montreal, montreal_gas_p[:,:rT])
	update_param!(m, :montreal_cycles, :rA_montreal, montreal_gas_p[:,:rA])
	update_param!(m, :montreal_cycles, :GU_montreal_0, montreal_init[:, :cumulative_uptake])
	update_param!(m, :montreal_cycles, :g0_montreal, g0_montreal)
	update_param!(m, :montreal_cycles, :g1_montreal, g1_montreal)
	update_param!(m, :montreal_cycles, :R0_montreal, Array(montreal_init[:,[:pool1,:pool2,:pool3,:pool4]]))
	update_param!(m, :montreal_cycles, :emiss2conc_montreal, montreal_gas_p[:,:emis2conc])
	update_param!(m, :montreal_cycles, :a_montreal, montreal_a)
	update_param!(m, :montreal_cycles, :τ_montreal, montreal_τ)

	# ---- Tropospheric Ozone Precursors, Aerosols, & Reactive Gas Cycles (Aerosol+) ---- #
	update_param!(m, :aerosol_plus_cycles, :aerosol_plus_0, aerosol_plus_init[:, :concentration])
	update_param!(m, :aerosol_plus_cycles, :r0_aerosol_plus, aerosol_plus_gas_p[:,:r0])
	update_param!(m, :aerosol_plus_cycles, :rU_aerosol_plus, aerosol_plus_gas_p[:,:rC])
	update_param!(m, :aerosol_plus_cycles, :rT_aerosol_plus, aerosol_plus_gas_p[:,:rT])
	update_param!(m, :aerosol_plus_cycles, :rA_aerosol_plus, aerosol_plus_gas_p[:,:rA])
	update_param!(m, :aerosol_plus_cycles, :GU_aerosol_plus_0, aerosol_plus_init[:, :cumulative_uptake])
	update_param!(m, :aerosol_plus_cycles, :g0_aerosol_plus, g0_aerosol_plus)
	update_param!(m, :aerosol_plus_cycles, :g1_aerosol_plus, g1_aerosol_plus)
	update_param!(m, :aerosol_plus_cycles, :R0_aerosol_plus, Array(aerosol_plus_init[:,[:pool1,:pool2,:pool3,:pool4]]))
	update_param!(m, :aerosol_plus_cycles, :emiss2conc_aerosol_plus, aerosol_plus_gas_p[:,:emis2conc])
	update_param!(m, :aerosol_plus_cycles, :a_aerosol_plus, aerosol_plus_a)
	update_param!(m, :aerosol_plus_cycles, :τ_aerosol_plus, aerosol_plus_τ)

	# ---- Radiative Forcing ---- #
	update_param!(m, :radiative_forcing, :co2_f, vec(Array(co2_p[:, [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :ch4_f, vec(Array(ch4_p[:, [:f1, :f2, :f3]])))
   	update_param!(m, :radiative_forcing, :ch4_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="methane|o3"), [:f1, :f2, :f3]])) )
   	update_param!(m, :radiative_forcing, :ch4_h2o_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="methane|strat_h2o"), [:f1, :f2, :f3]])) )
  	update_param!(m, :radiative_forcing, :n2o_f, vec(Array(n2o_p[:, [:f1, :f2, :f3]])))
   	update_param!(m, :radiative_forcing, :n2o_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nitrous_oxide|o3"), [:f1, :f2, :f3]])) )
	update_param!(m, :radiative_forcing, :montreal_f, Array(montreal_gas_p[:,[:f1, :f2, :f3]]))
	update_param!(m, :radiative_forcing, :montreal_ind_f, Array(montreal_irf_p[:,[:f1, :f2, :f3]]))
	update_param!(m, :radiative_forcing, :flourinated_f, Array(flourinated_gas_p[:,[:f1, :f2, :f3]]))
	update_param!(m, :radiative_forcing, :bc_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="bc"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :bc_snow_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="bc|bc_on_snow"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :bc_aci_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="bc|aci"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :co_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="co"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :co_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="co|o3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nh3_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nh3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nmvoc_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nmvoc"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nmvoc_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nmvoc|o3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nox"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_o3_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nox|o3"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_avi_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="nox_avi"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :nox_avi_contrails_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="nox_avi|contrails"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :oc_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="oc"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :oc_aci_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="oc|aci"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :so2_f, vec(Array(aerosol_plus_gas_p[(aerosol_plus_gas_p.gas_name.=="so2"), [:f1, :f2, :f3]])))
	update_param!(m, :radiative_forcing, :so2_aci_f, vec(Array(irf_p[(irf_p.indirect_forcing_effect.=="so2|aci"), [:f1, :f2, :f3]])))
  	update_param!(m, :radiative_forcing, :solar_f, 1.0)
  	update_param!(m, :radiative_forcing, :landuse_f, 1.0)
  	update_param!(m, :radiative_forcing, :volcanic_f, 1.0)
  	update_param!(m, :radiative_forcing, :solar_RF, forcing_data[:, :solar])
  	update_param!(m, :radiative_forcing, :landuse_RF, forcing_data[:, :land_use])
  	update_param!(m, :radiative_forcing, :volcanic_RF, forcing_data[:, :volcanic])
  	update_param!(m, :radiative_forcing, :other_RF, zeros(length(start_year:end_year)))

    # ---- Global Temperature Anomaly ---- #
    update_param!(m, :temperature, :Tj_0, vec(Array(init_thermal_vals[:,[:thermal_box_1,:thermal_box_2,:thermal_box_3]])))
    update_param!(m, :temperature, :T_0, init_thermal_vals.global_temp_anomaly[1])
    update_param!(m, :temperature, :q, vec(Array(thermal_p[:,[:q1,:q2,:q3]])))
    update_param!(m, :temperature, :decay_factor, thermal_decay_factors)

    # --- Set Parameters Common to Multiple Components ---- #
	add_shared_param!(m, :shared_co2_pi, co2_p.PI_conc[1])
    connect_param!(m, :co2_cycle, :co2_pi, :shared_co2_pi)
    connect_param!(m, :radiative_forcing, :co2_pi, :shared_co2_pi)

	add_shared_param!(m, :shared_ch4_pi, ch4_p.PI_conc[1])
    connect_param!(m, :ch4_cycle, :ch4_pi, :shared_ch4_pi)
    connect_param!(m, :radiative_forcing, :ch4_pi, :shared_ch4_pi)

	add_shared_param!(m, :shared_n2o_pi, n2o_p.PI_conc[1])
    connect_param!(m, :n2o_cycle, :n2o_pi, :shared_n2o_pi)
    connect_param!(m, :radiative_forcing, :n2o_pi, :shared_n2o_pi)

	add_shared_param!(m, :shared_montreal_pi, montreal_gas_p[:,:PI_conc], dims = [:montreal_gases])
    connect_param!(m, :montreal_cycles, :montreal_pi, :shared_montreal_pi)
    connect_param!(m, :radiative_forcing, :montreal_pi, :shared_montreal_pi)

	add_shared_param!(m, :shared_flourinated_pi, flourinated_gas_p[:,:PI_conc], dims = [:flourinated_gases])
    connect_param!(m, :flourinated_cycles, :flourinated_pi, :shared_flourinated_pi)
    connect_param!(m, :radiative_forcing, :flourinated_pi, :shared_flourinated_pi)

	add_shared_param!(m, :shared_aerosol_plus_pi, aerosol_plus_gas_p[:,:PI_conc], dims = [:aerosol_plus_gases])
    connect_param!(m, :aerosol_plus_cycles, :aerosol_plus_pi, :shared_aerosol_plus_pi)
    connect_param!(m, :radiative_forcing, :aerosol_plus_pi, :shared_aerosol_plus_pi)

 	# ---------------------------------------------
    # Create Connections Between Mimi Components
    # ---------------------------------------------

    # Syntax is :component_needing_a_parameter_input => :name_of_that_parameter, :component_calculating_required_values => :name_of_variable_output
    connect_param!(m, :co2_cycle           => :E_co2,             :emissions           => :E_co2)
    connect_param!(m, :ch4_cycle           => :E_ch4,             :emissions           => :E_ch4)
    connect_param!(m, :n2o_cycle           => :E_n2o,             :emissions           => :E_n2o)
    connect_param!(m, :flourinated_cycles  => :E_flourinated,     :emissions           => :E_flourinated)
    connect_param!(m, :montreal_cycles     => :E_montreal,        :emissions           => :E_montreal)
    connect_param!(m, :aerosol_plus_cycles => :E_aerosol_plus,    :emissions           => :E_aerosol_plus)
    connect_param!(m, :montreal_cycles     => :Tj,                :temperature         => :Tj)
    connect_param!(m, :flourinated_cycles  => :Tj,                :temperature         => :Tj)
    connect_param!(m, :aerosol_plus_cycles => :Tj,                :temperature         => :Tj)
    connect_param!(m, :co2_cycle           => :Tj,                :temperature         => :Tj)
    connect_param!(m, :ch4_cycle           => :Tj,                :temperature         => :Tj)
    connect_param!(m, :n2o_cycle           => :Tj,                :temperature         => :Tj)
    connect_param!(m, :radiative_forcing   => :co2_conc,          :co2_cycle     	   => :co2)
    connect_param!(m, :radiative_forcing   => :ch4_conc,          :ch4_cycle     	   => :ch4)
    connect_param!(m, :radiative_forcing   => :n2o_conc,          :n2o_cycle     	   => :n2o)
    connect_param!(m, :radiative_forcing   => :montreal_conc,     :montreal_cycles     => :montreal_conc)
    connect_param!(m, :radiative_forcing   => :flourinated_conc,  :flourinated_cycles  => :flourinated_conc)
    connect_param!(m, :radiative_forcing   => :aerosol_plus_conc, :aerosol_plus_cycles => :aerosol_plus_conc)
    connect_param!(m, :temperature   	   => :F,                 :radiative_forcing   => :total_RF)
end

