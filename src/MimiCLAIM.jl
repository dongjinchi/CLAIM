using Mimi

include("components_climate/emissions.jl")
include("components_climate/co2_cycle.jl")
include("components_climate/ch4_cycle.jl")
include("components_climate/n2o_cycle.jl")
include("components_climate/montreal_gas_cycles.jl")
include("components_climate/flourinated_gas_cycles.jl")
include("components_climate/aerosol_plus_gas_cycles.jl")
include("components_climate/radiative_forcing.jl")
include("components_climate/temperature.jl")
include("components_socio/SocioEconomicComponent.jl")
include("components_climate/ClimateRegionalComponent.jl")
include("components_damage/claim_damage.jl")
include("components_damage/bhm_damage.jl")
include("helper/helper_global.jl")
include("helper/helper_damage.jl")
include("helper/helper_fairv2.jl")
include("helper/helper_scc.jl")
include("helper/helper_mesmer.jl")
include("helper/helper_mc.jl")

# ============================================================================
# 【get_model】国家层级 B&K CLAIM 模型主入口。
# 功能：创建 Mimi 模型，设置 time/regions/horizon_dim，装配气候链 +
#       socioeconomic、climateregional、claim_damage 三组件。
# helper 分工：helper_global.jl 负责全局层（维度/装配/连线/shared 参数/
#       socio 与 climateregional 参数）；helper_damage.jl 只负责 claim_damage
#       损害组件参数；helper_fairv2.jl 负责 FaIR 气候链参数；
#       helper_scc.jl 提供 compute_scco2（CO₂ 脉冲 + Ramsey 折现）。
# 气候模式（climate 关键字）：
#   :fair （默认）装配 FaIR v2.0 九组件链（emissions → 6 气体循环 → RF →
#          temperature），全球温度距平由 RCMIP `scenario` 排放情景内生产生，
#          temperature:T 连接 climateregional.inputtemp 与 claim_damage.T_global；
#          T0_global 传 NaN，由 claim_damage 在运行时取 damage_start_year 当年
#          温度（SCC 计算时由 compute_scco2 覆写为 Base/Pulse 共用参考温度）。
#          初始条件文件为 1950 年口径，故默认 start_year=1950 完成 spin-up。
#   :csv  沿用 shared_parameters/climate_T_global.csv 的外生温度路径
#         （历史值平推，仅用于损害模块的独立试运行与 smoke 测试）。
# 情景（scenario 关键字）：统一选择排放路径与社会经济基线，仅限
#         data/fair_data 实际具备的五种 RCMIP 情景：
#         ssp119/ssp126/ssp245/ssp370/ssp585。首位数字即 SSP 社会经济情景
#         （ssp119、ssp126→SSP1，ssp245→SSP2，ssp370→SSP3，ssp585→SSP5），
#         据此加载 IIASA SSP v3.2 的 population_sspX.csv /
#         socioeconomic_ypcgrowth_sspX.csv（PWT 11.0 历史 1950–2024 +
#         情景 2025–2100 + DICE 式外推 2101–2500）。climate=:csv 模式下
#         排放路径不生效，但社会经济基线仍按 scenario 选择。
# 时间轴：默认 end_year=2150；RCMIP 情景数据到 2500 年，因此 :fair 模式下
#         end_year 可取 2500 以前的任何一年。
# 损害起始：默认 bau_ref_year=2025（FaIR 仍从 1950 spin-up；SCC 脉冲年默认
#         与此对齐为 2025，积分至 end_year）。
# 交互式 SCC / 出图：见 03Claim/MimiCLAIM.ipynb（notebook 内含本文件全部入口代码）。
# ============================================================================
function get_model(resolution::Symbol; climate::Symbol=:fair, scenario::String="ssp245",
                   damage_function::Symbol=:BK_region, income_state::Symbol=:dynamic,
                   capital_response::Bool=true,
                   theta_id::Int=0, start_year::Int=(climate == :fair ? 1950 : 2020),
                   end_year::Int=2150,
                   data_dir::AbstractString=joinpath(@__DIR__, "..", "data"), bau_ref_year::Int=2025,
                   TCR::Float64=1.79, RWF::Float64=0.552, F2x::Float64=3.759,
                   allow_legacy_climate_anchor::Bool=false)
    resolution == :national || error("Only resolution=:national is supported by MimiCLAIM.jl")
    climate in (:fair, :csv) || error("climate must be :fair or :csv")
    scenario in ("ssp119", "ssp126", "ssp245", "ssp370", "ssp585") ||
        error("scenario must be one of ssp119/ssp126/ssp245/ssp370/ssp585")
    validate_damage_controls(damage_function, income_state, capital_response)
    damage_component = active_damage_component(damage_function)
    ssp = _socio_ssp_from_scenario(scenario)
    start_year < end_year || error("start_year must be earlier than end_year")
    start_year <= bau_ref_year <= end_year || error("bau_ref_year must lie inside the model horizon")
    if climate == :fair
        start_year == 1950 || error("climate=:fair currently requires start_year=1950 because FaIR initial conditions are calibrated to 1950")
        end_year <= 2500 || error("climate=:fair requires end_year <= 2500 (RCMIP data range)")
    end

    shared_dir = joinpath(data_dir, "shared_parameters")
    regions = load_national_regions(shared_dir)
    years = collect(start_year:end_year)

    m = Model()
    setup_claim_dimensions!(m, years, regions)
    if climate == :fair
        # R25-c：FaIR 原始输入随 data_dir 走，不再硬编码到源码树。
        set_fair_gas_dimensions!(m; fair_data_dir=joinpath(data_dir, "fair_data"))
        add_fair_components!(m)
    end
    add_claim_components!(m; damage_component=damage_component)
    add_claim_shared_params!(m; data_dir=data_dir, years=years, regions=regions, ssp=ssp)
    connect_claim_components!(m; damage_component=damage_component)
    if climate == :fair
        update_MimiFAIR2_params!(m; emissions_forcing_scenario=scenario,
                                 start_year=start_year, end_year=end_year,
                                 TCR=TCR, RWF=RWF, F2x=F2x,
                                 fair_data_dir=joinpath(data_dir, "fair_data"))
        connect_param!(m, :climateregional, :inputtemp, :temperature, :T)
        if damage_component == :claim_damage
            connect_param!(m, :claim_damage, :T_global, :temperature, :T)
        end
    end
    update_global_params!(m; data_dir=data_dir, years=years, regions=regions, climate=climate, ssp=ssp,
                          start_year=start_year,
                          allow_legacy_climate_anchor=allow_legacy_climate_anchor)
    if damage_component == :claim_damage
        update_claim_damage_params!(m; data_dir=data_dir, years=years, regions=regions,
                                  theta_id=theta_id, start_year=start_year, bau_ref_year=bau_ref_year,
                                  climate=climate, damage_function=damage_function,
                                  income_state=income_state, capital_response=capital_response)
    else
        update_bhm_damage_params!(m; data_dir=data_dir, years=years, regions=regions,
                                  start_year=start_year, bau_ref_year=bau_ref_year)
    end
    return m
end
