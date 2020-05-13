#yearly aggregates
result_simulation_yearly_aggregate = DataFrame(Annual_generaion_tWh = Float64[], Annual_cost_billion_INR = Float64[], Annual_Emissions_CO2_billion_tonnes = Float64[],
                                        Annual_Emissions_SO2_million_tonnes = Float64[], Unit_cost_INRperkWh = Float64[], Unit_CO2_emissions_gm_per_kWh = Float64[],
                                        Unit_SO2_emissions_gm_per_kWh = Float64[], Renewable_inc_hydro_GWh = Float64[], Renewable_exc_hydro_GWh = Float64[],
                                        Solar_Gwh = Float64[], Wind_GWh = Float64[], Starts = Int[], Start_cost_mil_INR = Float64[], Start_cost_percentage = Float64[], Spilling_GWh = Float64[], Non_supply = Float64[])
push!(result_simulation_aggregate,(sum(result_simulation_df[:,:Generation_tWh]), sum(result_simulation_df[:,:Total_cost_billion_INR]), sum(result_simulation_df[:,:Emissions_CO2_million_tonnes])/10^3,
      sum(result_simulation_df[:,:Emissions_SO2_thousand_tonnes])/10^3, mean(result_simulation_df[:,:Unit_cost_INRperkWh]),mean(result_simulation_df[:,:Unit_CO2_emissions_gm_per_kWh]),
      mean(result_simulation_df[:,:Unit_SO2_emissions_gm_per_kWh]),sum(result_simulation_df[:,:Renewable_inc_hydro_GWh]), sum(result_simulation_df[:,:Renewable_exc_hydro_GWh]),
      sum(result_simulation_df[:,:Solar_Gwh]), sum(result_simulation_df[:,:Wind_GWh]), sum(result_simulation_df[:,:Starts]),sum(result_simulation_df[:,:Start_cost_mil_INR]),
      sum(result_simulation_df[:,:Start_cost_mil_INR])/(10^3*sum(result_simulation_df[:,:Total_cost_billion_INR])), sum(result_simulation_df[:,:Spilling_MWh])/10^3,
      sum(result_simulation_df[:,:Non_supply_MWh])/10^3 ))

XLSX.writetable("result_simulation_df.xlsx", collect(DataFrames.eachcol(result_simulation_df)), DataFrames.names(result_simulation_df))
XLSX.writetable("result_simulation_cap_uti_df.xlsx", collect(DataFrames.eachcol(result_simulation_cap_uti_df)), DataFrames.names(result_simulation_cap_uti_df))
XLSX.writetable("result_simulation_aggregate.xlsx", collect(DataFrames.eachcol(result_simulation_yearly_aggregate)), DataFrames.names(result_simulation_yearly_aggregate))
