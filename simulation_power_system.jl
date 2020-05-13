# simulation loop
for day = 1:36

param_dem_dict = Dict()
for i in 1:size(dem,2)
    k = 1
    for j in (size(temporal_resolution,1)*(day-1)+1):size(temporal_resolution,1)*(day-1)+size(temporal_resolution,1)
       param_dem_dict[dem[1,i],k] = dem[j+1,i]
       k += 1
     end
   end

param_capacity_factor_dict = Dict()

let k =1
for i in capacity_factor[1:end,1]
  d = 1
   for j in  size(temporal_resolution,1)*(day-1)+1:size(temporal_resolution,1)*(day-1)+size(temporal_resolution,1)
    param_capacity_factor_dict[i,d] = capacity_factor[k,j+2]
    d += 1
  end
k += 1
end
end

power_system = Model(with_optimizer(Gurobi.Optimizer, TimeLimit = 3600))

@variable(power_system,output[r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution;
         haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i])  && haskey(param_fuel_gen_dict,[f,i])] >= 0)
@variable(power_system,reserve_var[r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution, p in reserve;
          haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i)] >= 0)
@variable(power_system, not_supply[j in demand_centers, k in temporal_resolution] == 0)
@variable(power_system, spilling[j in demand_centers, k in temporal_resolution] == 0)
@variable(power_system,sto[r in regions, i in storage, j in demand_centers, k in temporal_resolution, l in storage_top;
                 haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_sto_dict,[r,i]) ] >= 0)
@variable(power_system,com_st[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)] >= 0,Int)
@variable(power_system,start_event[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)] >= 0,Int)
@variable(power_system,start_event_type[i in generators, k in temporal_resolution, j in start; haskey(param_firm_cap_dict,i)] >= 0,Int)


@variable(power_system,shut_event[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)] >= 0,Int)

@constraint(power_system,demand_supply[j in demand_centers,k in temporal_resolution],
            sum(output[r,f,i,j,k] * (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) for i in generators, r in regions, f in fuels
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])) +
            sum(sto[r,s,j,k,n]*((n == "in") ? -1 : param_sto_eff_dict[s]) * (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) for r in regions, s in storage, n in storage_top
            if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r])) + not_supply[j,k] - spilling[j,k]
            == param_dem_dict[j,k])
@constraint(power_system,reserve_const[j in demand_centers, k in temporal_resolution, p in reserve ],
            sum(reserve_var[r,f,i,j,k,p] for i in generators, r in regions, f in fuels
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i))
             >= .03 * param_dem_dict[j,k])

@constraint(power_system,capacity_supply[i in generators, k in temporal_resolution],
            sum(output[r,f,i,j,k] for r in regions , j in demand_centers, f in fuels
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            + sum(reserve_var[r,f,i,j,k,"up"] for r in regions, j in demand_centers, f in fuels
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i))
            <= (haskey(param_firm_cap_dict,i) ? 0.95 * param_cap_dict[i]*com_st[i,k] : param_cap_dict[i]*param_capacity_factor_dict[i,k]))

@constraint(power_system,sto_capacity[s in storage, k in temporal_resolution],
            sum(sto[r,s,j,p,n] * (n == "in" ? (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) : -1)  for r in regions, j in demand_centers, n in storage_top, p in 1:k
            if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r])) <= param_sto_cap_dict[s])
@constraint(power_system,sto_dispatch[s in storage, k in temporal_resolution; k !=1],
            sum(sto[r,s,j,k,"out"] for r in regions, j in demand_centers if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r]) ) <=
            sum(sto[r,s,j,p,n] * (n == "in" ? (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) : -1)  for r in regions, j in demand_centers, n in storage_top, p in 1:k-1
            if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r])))
@constraint(power_system,sto_dispatch_boundary[s in storage, k = 1],
            sum(sto[r,s,j,k,"out"] for r in regions, j in demand_centers if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r]) ) <=
            sum(sto[r,s,j,p,n] * (n == "in" ? (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) : -1)  for r in regions, j in demand_centers, n in storage_top, p = k
            if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r])))
@constraint(power_system,sto_rate[s in storage, k in temporal_resolution, n in storage_top],
            sum(sto[r,s,j,k,n] * (n == "in" ? (haskey(param_reg_dem_dict,[r,j]) ? 1/.95 : 1/.9) : 1)  for r in regions, j in demand_centers
                            if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r])) <= param_sto_dis_rate_dict[s])
@constraint(power_system,cluster_units[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
            com_st[i,k] <= param_gen_cluster_dict[i])
@constraint(power_system,stable_load[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
            sum(output[r,f,i,j,k] for r in regions , j in demand_centers, f in fuels
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                >= param_stable_load_dict[i] * com_st[i,k] + sum(reserve_var[r,f,i,j,k,"down"] for r in regions, j in demand_centers, f in fuels
                if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])))

@constraint(power_system,start_event_cons[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i) && k != 1],
            start_event[i,k] - shut_event[i,k] == com_st[i,k] - com_st[i,k-1])
@constraint(power_system,start_event_type_const[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i) && k >= 1],
                        start_event[i,k] == sum(start_event_type[i,k,j] for j in start))
@constraint(power_system, up_time_const[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
                        com_st[i,k] >= start_event[i,k])

if day == 1

@constraint(power_system,start_event_cons_boundary[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i) && k == 1],
            start_event[i,k] == com_st[i,k])
@constraint(power_system, down_time_cons[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
            param_gen_cluster_dict[i] - com_st[i,k] >= sum(shut_event[i,p] for p in (((k-param_gen_downtime_dict[i])>0) ? (k-param_gen_downtime_dict[i]) : 1) : k))
@constraint(power_system,start_event_beg_cons[i in generators, k in temporal_resolution, j in start; haskey(param_firm_cap_dict,i) && k ==1 && j != "cold"],
            start_event_type[i,k,j] == 0)
@constraint(power_system, start_event_hot_cons[i in generators, k in temporal_resolution, j in start; haskey(param_firm_cap_dict,i) && j == "hot"],
            start_event_type[i,k,j] <= (k < 12 ? sum(shut_event[i,l] for l in 1 : k) : sum(shut_event[i,l] for l in k-11 : k)))
@constraint(power_system, start_event_warm_cons[i in generators, k in temporal_resolution, j in start; haskey(param_firm_cap_dict,i) && j == "warm"],
              start_event_type[i,k,j] <= (k < 12 ? 0 : k < 16 ? (sum(shut_event[i,l] for l in 1 : k-11)) : (sum(shut_event[i,l] for l in k-15 : k - 12))))

else

@constraint(power_system,start_event_cons_boundary[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i) && k == 1],
            start_event[i,k] - shut_event[i,k] == com_st[i,k] - result_com_st_dict[[i,size(temporal_resolution,1)]])
@constraint(power_system, down_time_cons[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
            param_gen_cluster_dict[i] - com_st[i,k] >=
            ((k-param_gen_downtime_dict[i] <= 0)
            ? sum(result_shut_event_lag_one_day_dict[[i,p]] for p in (size(temporal_resolution,1) - param_gen_downtime_dict[i] + k) : size(temporal_resolution,1)) +
            sum(shut_event[i,p] for p in 1 : k)
            : sum(shut_event[i,p] for p in (k-param_gen_downtime_dict[i]) : k)))
@constraint(power_system, start_event_hot_cons[i in generators, k in temporal_resolution, j in start; haskey(param_firm_cap_dict,i) && j == "hot" && k != size(temporal_resolution,1)],
            start_event_type[i,k,j] <= (k < 12 ? (sum(shut_event[i,l] for l in 1 : k) + sum(result_shut_event_lag_one_day_dict[[i,l]] for l in (size(temporal_resolution,1)-11+k)
            : size(temporal_resolution,1))) : sum(shut_event[i,l] for l in k-11 : k)))
@constraint(power_system,start_event_warm_cons[i in generators, k in temporal_resolution, j in start; haskey(param_firm_cap_dict,i) && j == "warm"],
            start_event_type[i,k,j] <= (k < 12 ? (sum(result_shut_event_lag_one_day_dict[[i,l]] for l in size(temporal_resolution,1)-16+k : size(temporal_resolution,1)-16+k+4)) :
            (k > 16) ? sum(shut_event[i,l] for l in k-15 : k-11) : (sum(result_shut_event_lag_one_day_dict[[i,l]] for l in size(temporal_resolution,1)-16+k : size(temporal_resolution,1)) + sum(shut_event[i,l] for l in 1 : k-11))
            ))

@constraint(power_system,ramp_down_roll[i in generators, k = 1 ; haskey(param_firm_cap_dict,i)],
            sum(result_output_dict[[r,f,i,j,size(temporal_resolution,1)]] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            -sum(output[r,f,i,j,k] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            <= result_com_st_dict[[i,size(temporal_resolution,1)]]*param_ramp_down_dict[i] + param_stable_load_dict[i]*shut_event[i,k])

@constraint(power_system,ramp_up_roll[i in generators, k = 1 ; haskey(param_firm_cap_dict,i)],
            sum(result_output_dict[[r,f,i,j,size(temporal_resolution,1)]] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            -sum(output[r,f,i,j,k] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            >= -result_com_st_dict[[i,size(temporal_resolution,1)]]*param_ramp_up_dict[i] - param_stable_load_dict[i]*start_event[i,k])
end

@constraint(power_system,ramp_down[i in generators,k in temporal_resolution; haskey(param_firm_cap_dict,i) && k != temporal_resolution[length(temporal_resolution)]],
            sum(output[r,f,i,j,k] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            -sum(output[r,f,i,j,k+1] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            <= com_st[i,k]*param_ramp_down_dict[i] + param_stable_load_dict[i]*shut_event[i,k+1])
@constraint(power_system,ramp_up[i in generators,k in temporal_resolution; haskey(param_firm_cap_dict,i) && k != temporal_resolution[length(temporal_resolution)]],
            sum(output[r,f,i,j,k] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            -sum(output[r,f,i,j,k+1] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
              >= -com_st[i,k]*param_ramp_up_dict[i] - param_stable_load_dict[i]*start_event[i,k+1])
@constraint(power_system,transmission_constraint[r in regions, p in regions, k in temporal_resolution; haskey(param_reg_exc_dict,[r,p])],
            sum(output[r,f,i,j,k] for i in generators, f in fuels, j in demand_centers
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_dem_dict,[p,j]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])) +
            sum(sto[r,s,j,k,n]*((n == "in") ? -1 : param_sto_eff_dict[s]) * (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) for r in regions, s in storage, n in storage_top, j in demand_centers
            if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_dem_dict,[p,j]))
            <= param_reg_exc_dict[[r,p]])

# Cost minimizing
@objective(power_system, Min, sum(output[r,f,i,j,k]*param_gen_cost_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
            if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                + sum(start_event_type[i,k,j]*param_startup_cost_dict[[i,j]] for r in regions, f in fuels, i in generators, k in temporal_resolution, j in start
                if haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i))
                + sum(sto[r,s,j,k,n]*param_sto_cost_dict[s] for r in regions, s in storage, n in ["out"], j in demand_centers, k in temporal_resolution
                    if haskey(param_reg_sto_dict,[r,s]) && haskey(param_dem_top_dict,[j,r]))
                + sum(not_supply[j,k]*100*param_non_supply_dict[j] for j in demand_centers, k in temporal_resolution)
                + sum(spilling[j,k]*1000 for j in demand_centers, k in temporal_resolution)
                + sum(reserve_var[r,f,i,j,k,p]*param_gen_cost_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution, p in reserve
                if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i)))

# Emission minimizing
#@objective(power_system, Min, sum(output[r,f,i,j,k]*param_gen_emissions_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
    #        if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])))

optimize!(power_system)

#print(power_system)
# writing results to excel by first converting it to a DataFrame


result_output_df = DataFrame(Region = String[], Fuel = String[], Generator = String[], Demand_center = String[], Time = Int[], Generation = Float64[])
result_demand_met_df = DataFrame(Region = String[], Time = Int[], Demand_met = Float64[])
result_reserve_df = DataFrame(Region = String[], Fuel = String[], Generator = String[], Demand_center = String[], Time = Int[], Direction = String[], Contribution = Float64[])
result_spilling_df = DataFrame(Region = String[], Time = Int[], spilled = Float64[])
result_non_supply_df = DataFrame(Region = String[], Time = Int[], non_supply = Float64[])
result_start_event_df = DataFrame(Generator = String[], Time = Int[], start_event = Int[])
result_shut_event_df = DataFrame(Generator = String[], Time = Int[], shut_event = Int[])
result_com_st_df = DataFrame(Generator = String[], Time = Int[], Committed = Int[])
result_output_transmission_df = DataFrame(From = String[], To = String[], Time = Int[], Generation = Float64[])

result_start_shut_com_df = DataFrame(Generator = String[], Time = Int[], Start = Float64[], Shut = Float64[], Com = Float64[], Cold = Float64[], Warm = Float64[], Hot = Float64[])

for j in demand_centers , k in temporal_resolution
    push!(result_spilling_df,(j,k,value(spilling[j,k])))
    push!(result_non_supply_df,(j,k,value(not_supply[j,k])))
    push!(result_demand_met_df,(j,k,value(demand_supply[j,k])))
end

for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])
        push!(result_output_df,(r,f,i,j,k,value(output[r,f,i,j,k])))
        result_output_dict[[r,f,i,j,k]] = value(output[r,f,i,j,k])
        #value(output[i,j,k])
    end
end

for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution, p in reserve
    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i)
        push!(result_reserve_df,(r,f,i,j,k,p,value(reserve_var[r,f,i,j,k,p])))
        #value(output[i,j,k])
    end
end


for i in generators, k in temporal_resolution
    if haskey(param_firm_cap_dict,i)
        #push!(result_start_event_df,(i,k,value(start_event[i,k])))
        #push!(result_com_st_df,(i,k,value(com_st[i,k])))
        #push!(result_shut_event_df,(i,k,value(shut_event[i,k])))
        push!(result_start_shut_com_df,
        (i,k,value(start_event[i,k]),value(shut_event[i,k]),value(com_st[i,k]), value(start_event_type[i,k,"cold"]),value(start_event_type[i,k,"warm"]),
        value(start_event_type[i,k,"hot"])))
    end
end

for r in regions, p in regions, k in temporal_resolution
    if haskey(param_reg_exc_dict,[r,p])
    push!(result_output_transmission_df,(r,p,k,sum(value(output[r,f,i,j,k]) for i in generators, j in demand_centers, f in fuels
     if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_reg_dem_dict,[p,j]))))
end
end

for i in generators, j in temporal_resolution
    if haskey(param_firm_cap_dict,i)
        if day > 1
            result_shut_event_lag_two_day_dict[[i,j]]= result_shut_event_lag_one_day_dict[[i,j]]
        end
        result_shut_event_lag_one_day_dict[[i,j]] = value(shut_event[i,j])
        result_start_event_dict[[i,j]] = value(start_event[i,j])
        result_com_st_dict[[i,j]] = value(com_st[i,j])

    end
end
#
if day > 1
XLSX.writetable(string("result_output_df" , "$day" , ".xlsx"), collect(DataFrames.eachcol(result_output_df)), DataFrames.names(result_output_df))
XLSX.writetable(string("result_reserve_df" , "$day" , ".xlsx"), collect(DataFrames.eachcol(result_reserve_df)), DataFrames.names(result_reserve_df))
XLSX.writetable(string("result_output_transmission_df" , "$day" , ".xlsx"), collect(DataFrames.eachcol(result_output_transmission_df)), DataFrames.names(result_output_transmission_df))
#XLSX.writetable("result_start_event_df.xlsx", collect(DataFrames.eachcol(result_start_event_df)), DataFrames.names(result_start_event_df))
XLSX.writetable(string("result_start_shut_com_df" , "$day" , ".xlsx"), collect(DataFrames.eachcol(result_start_shut_com_df)), DataFrames.names(result_start_shut_com_df))
XLSX.writetable(string("result_spilling_df" , "$day", ".xlsx"), collect(DataFrames.eachcol(result_spilling_df)), DataFrames.names(result_spilling_df))
XLSX.writetable(string("result_non_supply_df" , "$day", ".xlsx"), collect(DataFrames.eachcol(result_non_supply_df)), DataFrames.names(result_non_supply_df))
XLSX.writetable(string("result_demand_met_df" , "$day", ".xlsx"), collect(DataFrames.eachcol(result_demand_met_df)), DataFrames.names(result_demand_met_df))
end

Generation = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
Emissions_CO2 = sum(value(output[r,f,i,j,k])*param_gen_emissions_CO2_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
Emissions_SO2 = sum(value(output[r,f,i,j,k])*param_gen_emissions_SO2_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
Renewable_inc_hydro = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                      if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Hydro" || f == "Wind" || f == "Solar"))
Renewable_exc_hydro = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Wind" || f == "Solar"))
Solar = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                     if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Solar"))
Wind = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Wind" ))
Coal = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Coal" ))
Gas = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Gas" ))
Nuclear = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Nuclear" ))
Diesel = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
                    if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && (f == "Diesel" ))


Starts = sum(value(start_event[i,k]) for i in generators, k in temporal_resolution if haskey(param_firm_cap_dict,i))
Start_cost = sum(value(start_event_type[i,k,j])*param_startup_cost_dict[[i,j]] for i in generators, k in temporal_resolution, j in start
            if haskey(param_firm_cap_dict,i))
Spilling = sum(value(spilling[j,k]) for j in demand_centers, k in temporal_resolution)
Non_supply = sum(value(not_supply[j,k]) for j in demand_centers, k in temporal_resolution)
Total_cost = objective_value(power_system)

push!(result_simulation_df,(day,Generation/10^6,Total_cost/10^9,Emissions_CO2/10^6, Emissions_SO2/10^9,Total_cost/(Generation*10^3),
      Emissions_CO2/Generation, Emissions_SO2/(Generation*10^3), Renewable_inc_hydro/10^3, Renewable_exc_hydro/10^3, Solar/10^3, Wind/10^3, Coal/10^3, Gas/10^3, Nuclear/10^3, Diesel/10^3,
       round(Starts), Start_cost/10^6, Start_cost/Total_cost, Spilling, Non_supply, MOI.get(power_system,MOI.RelativeGap()), MOI.get(power_system, MOI.RawStatusString())))


cap_uti = zeros(size(generators,1))
let p =1
for i in generators
    cap_uti[p] = sum(value(output[r,f,i,j,k]) for r in regions, f in fuels, j in demand_centers, k in temporal_resolution
                  if haskey(param_dem_top_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))/(param_cap_dict[i]*param_gen_cluster_dict[i]*size(temporal_resolution,1))
    p +=1
end
end


insertcols!(result_simulation_cap_uti_df,day+1,name = cap_uti, makeunique = true)

end
