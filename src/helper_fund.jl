using Distributions

"""
Converts a year value into an integer corresponding to fund's time index.
"""
function getindexfromyear(year::Int)
    baseyear = 1990
    return year - baseyear + 1
end

"""
Reads parameter csvs from data directory into a dictionary with two keys:
* :shared => (parameter_name => default_value) for parameters shared in the model
* :unshared => ((component_name, parameter_name) => default_value) for component specific parameters that are not shared

For parameters defined as distributions, this sets the value to their mode.
""" 
function load_default_parameters(nsteps)
    datadir = joinpath(dirname(@__FILE__), "..", "data")
    # Load the unshared parameters
    files = readdir(joinpath(datadir, "unshared_parameters"))
    filter!(i -> (i != "desktop.ini" && i != ".DS_Store"), files)
    unshared_parameters = Dict{Tuple{Symbol,Symbol}, Any}()
    for file in files
        param_info = Symbol.(split(splitext(file)[1], "-"))
        unshared_parameters[(param_info[1], param_info[2])] = readdlm(joinpath(datadir,"unshared_parameters",file), ',', comments = true)
    end
    prepparameters!(unshared_parameters, nsteps)

    # Handle the shared parameters
    files = readdir(joinpath(datadir, "shared_parameters"))
    filter!(i -> (i != "desktop.ini" && i != ".DS_Store"), files)
    shared_parameters = Dict{Symbol, Any}(Symbol(splitext(file)[1]) => readdlm(joinpath(datadir,"shared_parameters",file), ',', comments = true) for file in files)
    prepparameters!(shared_parameters, nsteps)

    return Dict(:shared => shared_parameters, :unshared => unshared_parameters)
end


# For Truncated Gamma distributions, fund uses the mode of the untrucated distribution as it's default value.
import StatsBase.mode
function mode(d::Truncated{Gamma{Float64},Continuous})
    return mode(d.untruncated)
end


"""
Returns the mode for a distributional parameter; returns the value if it's not a distribution.
"""
getbestguess(p) = isa(p, ContinuousUnivariateDistribution) ? mode(p) : p


"""
Converts the original parameter dictionary loaded from the data files into a dictionary of default parameter values.

Original dictionary: parameter_name => string of distributions or values from csv file
Final dictionary: parameter_name => default value
"""
function prepparameters!(parameters, nsteps)
    for (param_name, p) in parameters
        column_count = size(p,2)
        if column_count == 1
            parameters[param_name] = getbestguess(convertparametervalue(p[1,1]))
        elseif column_count == 2
            #lengthp = size(p,1) #<= nsteps ? size(p,1) : nsteps
            parameters[param_name] = Float64[getbestguess(convertparametervalue(p[j,2])) for j in 1:size(p,1)]
        elseif column_count == 3
            length_index1 = length(unique(p[:,1])) <= nsteps ? length(unique(p[:,1])) : nsteps
            length_index2 = length(unique(p[:,2]))
            new_p = Array{Float64}(undef, length_index1, length_index2)
            cur_1 = 1
            cur_2 = 1
            lengthp = length_index1 * length_index2
            for j in 1:lengthp #size(p,1) 
                new_p[cur_1,cur_2] = getbestguess(convertparametervalue(p[j,3]))
                cur_2 += 1
                if cur_2 > length_index2
                    cur_2 = 1
                    cur_1 += 1
                end
            end
            parameters[param_name] = new_p
        end
    end
end


"""
Takes as input a single parameter value. 
If the parameter value is a string containing a distribution definition, it returns the distribtion.
If the parameter value is a number, it returns the number.
"""
function convertparametervalue(pv)
    if isa(pv,AbstractString)
        if startswith(pv,"~") && endswith(pv,")")
            args_start_index = something(findfirst(isequal('('), pv), 0) 
            dist_name = pv[2:args_start_index-1]
            args = split(pv[args_start_index+1:end-1], ';')
            fixedargs = filter(i->!occursin("=", i),args)
            optargs = Dict(split(i,'=')[1]=>split(i,'=')[2] for i in filter(i->occursin("=", i),args))

            if dist_name == "N"
                if length(fixedargs)!=2 error() end
                if length(optargs)>2 error() end

                basenormal = Normal(parse(Float64, fixedargs[1]),parse(Float64, fixedargs[2]))

                if length(optargs)==0
                    return basenormal
                else
                    return Truncated(basenormal,
                        haskey(optargs,"min") ? parse(Float64, optargs["min"]) : -Inf,
                        haskey(optargs,"max") ? parse(Float64, optargs["max"]) : Inf)
                end
            elseif startswith(pv, "~Gamma(")
                if length(fixedargs)!=2 error() end
                if length(optargs)>2 error() end

                basegamma = Gamma(parse(Float64, fixedargs[1]),parse(Float64, fixedargs[2]))

                if length(optargs)==0
                    return basegamma
                else
                    return Truncated(basegamma,
                        haskey(optargs,"min") ? parse(Float64, optargs["min"]) : -Inf,
                        haskey(optargs,"max") ? parse(Float64, optargs["max"]) : Inf)
                end
            elseif startswith(pv, "~Triangular(")
                triang = TriangularDist(parse(Float64, fixedargs[1]), parse(Float64, fixedargs[2]), parse(Float64, fixedargs[3]))
                return triang
            else
                error("Unknown distribution")
            end
        elseif pv=="true"
            return true
        elseif pv=="false"
            return false
        elseif endswith(pv, "y")
            return parse(Int, strip(pv,'y'))
        else
            try
                return parse(Float64, pv)
            catch e
                error(pv)
            end
        end
        return pv
    else
        return pv
    end
end

function update_fund_params!(m; nsteps)
    
    # ---------------------------------------------
    # Set component parameters
    # ---------------------------------------------
    parameters = load_default_parameters(nsteps)
    
    # Set unshared parameters - name is a Tuple{Symbol, Symbol} of (component_name, param_name)
    for (name, value) in parameters[:unshared]
        update_param!(m, name[1], name[2], value)
    end
        
    update_param!(m, :impactaggregation,:linearterm,  0.0162)
    update_param!(m, :impactaggregation,:quadrterm,  -0.000903)
        
    # ---------------------------------------------
    # Create Connections Between Mimi Components
    # ---------------------------------------------
    connect_param!(m, :impactaggregation, :regtmp,:climateregional, :regtmp)
    connect_param!(m, :impactaggregation, :temp2019,:climateregional, :temp2019)
    connect_param!(m, :impactaggregation, :income,:socioeconomic, :income)
    connect_param!(m, :impactaggregation, :ypc,:socioeconomic, :ypc)
    
    add_shared_param!(m, :ypcgrowth, parameters[:shared][:ypcgrowth], dims=[:time, :regions])
    connect_param!(m, :socioeconomic, :ypcgrowth, :ypcgrowth)
    connect_param!(m, :impactaggregation, :ypcgrowth, :ypcgrowth)
    
end
    