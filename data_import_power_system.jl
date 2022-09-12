# the new one 
import XLSX
data=XLSX.readxlsx("data_India_power_system_COVID.xlsx")
# demand_centers = readxlsheet("data_India_power_system_COVID.xlsx","demand_centers")
# region_dem_top = readxlsheet("data_India_power_system_COVID.xlsx","region_dem_top")
demand_centers=data["demand_centers"][:]
region_dem_top=data["region_dem_top"][:]
# temporal_resolution = [1 : 1 : readxlsheet("data_India_power_system_COVID.xlsx","temporal_resolution");]
# temporal_resolution=[1:1: data["temporal_resolution"][1]]
temporal_resolution = [i for i in range(1,25)]
temporal_resolution = [1:1:24;]
# storage_data = readxlsheet("data_India_power_system_COVID.xlsx","storage")
storage_data=data["storage"][:]
storage = [string(storage_data[i,2]) for i in 2:size(storage_data,1)]
storage_top = ["in","out"]
param_reg_sto_dict = Dict([storage_data[i+1,1],storage_data[i+1,2]] => 1 for i in 1:size(storage_data,1)-1)
param_sto_cap_dict = Dict(storage_data[i+1,2] => storage_data[i+1,4] for i in 1:size(storage_data,1)-1)
param_sto_eff_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,6] for i in 1:size(storage_data,1)-1)
param_sto_cost_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,7] for i in 1:size(storage_data,1)-1)
param_sto_ramp_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,8] for i in 1:size(storage_data,1)-1)
param_sto_dis_rate_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,5] for i in 1:size(storage_data,1)-1)


# gen_data = readxlsheet("data_India_power_system_COVID.xlsx","gen_data")
gen_data=data["gen_data"][:]
generators = [string(gen_data[i,2], gen_data[i,4]) for i in 2:size(gen_data,1)]
regions = unique(gen_data[2:end,1])
fuels = unique(gen_data[2:end,3])
ramp = gen_data[2:end,8]
stable_load = gen_data[2:end,7]
start = ["hot","warm","cold"]
reserve = ["up","down"]
startup_cost = gen_data[2:end,6]
down_time = gen_data[2:end,9]
clusters = gen_data[2:end,12]


# reg_exc = readxlsheet("data_India_power_system_COVID.xlsx","reg_exc")
reg_exc=data["reg_exc"][:]
# dem = readxlsheet("data_India_power_system_COVID.xlsx","demand_data")
dem =data["demand_data"][:]
# non_supply = readxlsheet("data_India_power_system_COVID.xlsx","non_supply")
non_supply =data["non_supply"][:]

# capacity_factor = readxlsheet("data_India_power_system_COVID.xlsx","capacity_factor")
capacity_factor=data["capacity_factor"][:]
capacity_factor[1:end,1] = [string(capacity_factor[i,1],1) for i in 1:size(capacity_factor,1)]

param_gen_cost_dict = Dict(generators[i] => gen_data[i+1,10] for i in 1:size(gen_data,1)-1)
param_gen_emissions_SO2_dict = Dict(generators[i] => gen_data[i+1,13] for i in 1:size(gen_data,1)-1)
param_gen_emissions_CO2_dict = Dict(generators[i] => gen_data[i+1,14] for i in 1:size(gen_data,1)-1)
param_reg_gen_dict = Dict([gen_data[i+1,1],generators[i]] => 1 for i in 1:size(generators,1))
param_non_supply_dict = Dict(non_supply[i,1] => non_supply[i,2] for i in 2:size(non_supply,1))

param_fuel_gen_dict = Dict([gen_data[i+1,3],generators[i]] => 1 for i in 1:size(generators,1))
param_dem_top_dict = Dict([region_dem_top[i,2],region_dem_top[i,3]] => 1 for i in 1:size(region_dem_top,1))
param_reg_dem_dict = Dict([region_dem_top[i,1],region_dem_top[i,2]] => 1 for i in 1:size(region_dem_top,1))
param_reg_exc_dict = Dict(reg_exc[i+1,1:2] => reg_exc[i+1,3] for i in 1:size(reg_exc,1)-1)

param_cap_dict = Dict(generators[i] => gen_data[i+1,5] for i in 1 : size(generators,1))
param_firm_cap_dict = Dict(generators[i] => 1 for i in 1 : size(generators,1) if gen_data[i+1,11] == "firm")
param_stable_load_dict = Dict(generators[i] => stable_load[i] for i in 1:size(stable_load,1))
param_startup_cost_dict = Dict([generators[i],j] => (j == "warm" ? startup_cost[i] : (j == "cold" ? startup_cost[i]*1.3 : startup_cost[i]*.7)) for i in 1:size(startup_cost,1), j in start)
param_gen_cluster_dict = Dict(generators[i] => clusters[i] for i in 1:size(clusters,1))
param_gen_downtime_dict = Dict(generators[i] => down_time[i] for i in 1:size(down_time,1))


param_ramp_up_dict = Dict(generators[i] => ramp[i] for i in 1:size(ramp,1))
param_ramp_down_dict = Dict(generators[i] => ramp[i] for i in 1:size(ramp,1))


# intializing the dataframe/containers to hold the results
result_shut_event_lag_one_day_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
result_shut_event_lag_two_day_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
result_shut_event_lag_three_day_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))

result_start_event_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
result_com_st_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
result_output_dict = Dict{Array{Any,1},Float64}([r,f,i,j,k] => 1 for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i])  && haskey(param_fuel_gen_dict,[f,i]))

result_simulation_df = DataFrame(Day = Int[], Generation_tWh = Float64[], Total_cost_billion_INR = Float64[], Emissions_CO2_million_tonnes = Float64[], Emissions_SO2_thousand_tonnes = Float64[], Unit_cost_INRperkWh = Float64[],
                                       Unit_CO2_emissions_gm_per_kWh = Float64[], Unit_SO2_emissions_gm_per_kWh = Float64[], Renewable_inc_hydro_GWh = Float64[], Renewable_exc_hydro_GWh = Float64[],
                                       Solar_Gwh = Float64[], Wind_GWh = Float64[], Coal_GWh = Float64[], Gas_GWh = Float64[], Nuclear_GWh = Float64[], Diesel_GWh = Float64[],
                                       Starts = Int[], Start_cost_mil_INR = Float64[], Start_cost_percentage = Float64[],
                                       Spilling_MWh = Float64[], Non_supply_MWh = Float64[], Relative_gap = Float64[], Termination_status = String[])
result_simulation_cap_uti_df = DataFrame(Generators = generators)


# the old one:

# demand_centers = readxlsheet("data_India_power_system_COVID.xlsx","demand_centers")
# region_dem_top = readxlsheet("data_India_power_system_COVID.xlsx","region_dem_top")
# # hello

# temporal_resolution = [1 : 1 : readxlsheet("data_India_power_system_COVID.xlsx","temporal_resolution");]

# storage_data = readxlsheet("data_India_power_system_COVID.xlsx","storage")
# storage = [string(storage_data[i,2]) for i in 2:size(storage_data,1)]
# storage_top = ["in","out"]
# param_reg_sto_dict = Dict([storage_data[i+1,1],storage_data[i+1,2]] => 1 for i in 1:size(storage_data,1)-1)
# param_sto_cap_dict = Dict(storage_data[i+1,2] => storage_data[i+1,4] for i in 1:size(storage_data,1)-1)
# param_sto_eff_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,6] for i in 1:size(storage_data,1)-1)
# param_sto_cost_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,7] for i in 1:size(storage_data,1)-1)
# param_sto_ramp_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,8] for i in 1:size(storage_data,1)-1)
# param_sto_dis_rate_dict =  Dict(storage_data[i+1,2] => storage_data[i+1,5] for i in 1:size(storage_data,1)-1)


# gen_data = readxlsheet("data_India_power_system_COVID.xlsx","gen_data")
# generators = [string(gen_data[i,2], gen_data[i,4]) for i in 2:size(gen_data,1)]
# regions = unique(gen_data[2:end,1])
# fuels = unique(gen_data[2:end,3])
# ramp = gen_data[2:end,8]
# stable_load = gen_data[2:end,7]
# start = ["hot","warm","cold"]
# reserve = ["up","down"]
# startup_cost = gen_data[2:end,6]
# down_time = gen_data[2:end,9]
# clusters = gen_data[2:end,12]


# reg_exc = readxlsheet("data_India_power_system_COVID.xlsx","reg_exc")
# dem = readxlsheet("data_India_power_system_COVID.xlsx","demand_data")
# non_supply = readxlsheet("data_India_power_system_COVID.xlsx","non_supply")

# capacity_factor = readxlsheet("data_India_power_system_COVID.xlsx","capacity_factor")
# capacity_factor[1:end,1] = [string(capacity_factor[i,1],1.0) for i in 1:size(capacity_factor,1)]

# param_gen_cost_dict = Dict(generators[i] => gen_data[i+1,10] for i in 1:size(gen_data,1)-1)
# param_gen_emissions_SO2_dict = Dict(generators[i] => gen_data[i+1,13] for i in 1:size(gen_data,1)-1)
# param_gen_emissions_CO2_dict = Dict(generators[i] => gen_data[i+1,14] for i in 1:size(gen_data,1)-1)
# param_reg_gen_dict = Dict([gen_data[i+1,1],generators[i]] => 1 for i in 1:size(generators,1))
# param_non_supply_dict = Dict(non_supply[i,1] => non_supply[i,2] for i in 2:size(non_supply,1))

# param_fuel_gen_dict = Dict([gen_data[i+1,3],generators[i]] => 1 for i in 1:size(generators,1))
# param_dem_top_dict = Dict([region_dem_top[i,2],region_dem_top[i,3]] => 1 for i in 1:size(region_dem_top,1))
# param_reg_dem_dict = Dict([region_dem_top[i,1],region_dem_top[i,2]] => 1 for i in 1:size(region_dem_top,1))
# param_reg_exc_dict = Dict(reg_exc[i+1,1:2] => reg_exc[i+1,3] for i in 1:size(reg_exc,1)-1)

# param_cap_dict = Dict(generators[i] => gen_data[i+1,5] for i in 1 : size(generators,1))
# param_firm_cap_dict = Dict(generators[i] => 1 for i in 1 : size(generators,1) if gen_data[i+1,11] == "firm")
# param_stable_load_dict = Dict(generators[i] => stable_load[i] for i in 1:size(stable_load,1))
# param_startup_cost_dict = Dict([generators[i],j] => (j == "warm" ? startup_cost[i] : (j == "cold" ? startup_cost[i]*1.3 : startup_cost[i]*.7)) for i in 1:size(startup_cost,1), j in start)
# param_gen_cluster_dict = Dict(generators[i] => clusters[i] for i in 1:size(clusters,1))
# param_gen_downtime_dict = Dict(generators[i] => down_time[i] for i in 1:size(down_time,1))


# param_ramp_up_dict = Dict(generators[i] => ramp[i] for i in 1:size(ramp,1))
# param_ramp_down_dict = Dict(generators[i] => ramp[i] for i in 1:size(ramp,1))


# # intializing the dataframe/containers to hold the results
# result_shut_event_lag_one_day_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
# result_shut_event_lag_two_day_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
# result_shut_event_lag_three_day_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))

# result_start_event_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
# result_com_st_dict = Dict{Array{Any,1},Float64}([i,j] => 1 for i in generators, j in temporal_resolution if haskey(param_firm_cap_dict,i))
# result_output_dict = Dict{Array{Any,1},Float64}([r,f,i,j,k] => 1 for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i])  && haskey(param_fuel_gen_dict,[f,i]))

# result_simulation_df = DataFrame(Day = Int[], Generation_tWh = Float64[], Total_cost_billion_INR = Float64[], Emissions_CO2_million_tonnes = Float64[], Emissions_SO2_thousand_tonnes = Float64[], Unit_cost_INRperkWh = Float64[],
#                                        Unit_CO2_emissions_gm_per_kWh = Float64[], Unit_SO2_emissions_gm_per_kWh = Float64[], Renewable_inc_hydro_GWh = Float64[], Renewable_exc_hydro_GWh = Float64[],
#                                        Solar_Gwh = Float64[], Wind_GWh = Float64[], Coal_GWh = Float64[], Gas_GWh = Float64[], Nuclear_GWh = Float64[], Diesel_GWh = Float64[],
#                                        Starts = Int[], Start_cost_mil_INR = Float64[], Start_cost_percentage = Float64[],
#                                        Spilling_MWh = Float64[], Non_supply_MWh = Float64[], Relative_gap = Float64[], Termination_status = String[])
# result_simulation_cap_uti_df = DataFrame(Generators = generators)
