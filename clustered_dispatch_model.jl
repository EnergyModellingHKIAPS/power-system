import Pkg; Pkg.add("DataFrames"); Pkg.add("XLSX"); Pkg.add("GLPK"); Pkg.add("Gurobi")
using DataFrames
using ExcelReaders
using GLPK
using XLSX
using Gurobi
using JuMP

demand_centers = readxlsheet("data_INDIA_clustered.xlsx","demand_centers")
region_dem_incoming = readxlsheet("data_INDIA_clustered.xlsx","region_dem_incoming")

#temporal_resolution = readxlsheet("data_INDIA_clustered.xlsx","temporal_resolution")

temporal_resolution = [1 : 1 : readxlsheet("data_INDIA_clustered.xlsx","temporal_resolution");]

gen_data = readxlsheet("data_INDIA_clustered.xlsx","gen_data")
generators = [string(gen_data[i,2], gen_data[i,4]) for i in 2:size(gen_data,1)]
regions = unique(gen_data[2:end,1])
fuels = unique(gen_data[2:end,3])
ramp = gen_data[2:end,8]
stable_load = gen_data[2:end,7]
startup_cost = gen_data[2:end,6]
clusters = gen_data[2:end,11]
reg_exc = readxlsheet("data_INDIA_clustered.xlsx","reg_exc")
dem = readxlsheet("data_INDIA_clustered.xlsx","demand_data")

capacity_factor = readxlsheet("data_INDIA_clustered.xlsx","capacity_factor")

param_gen_cost_dict = Dict(generators[i] => gen_data[i+1,9] for i in 1:size(gen_data,1)-1)
param_reg_gen_dict = Dict([gen_data[i+1,1],generators[i]] => 1 for i in 1:size(generators,1))
param_fuel_gen_dict = Dict([gen_data[i+1,3],generators[i]] => 1 for i in 1:size(generators,1))
param_dem_incoming_dict = Dict([region_dem_incoming[i,2],region_dem_incoming[i,3]] => 1 for i in 1:size(region_dem_incoming,1))
param_reg_dem_dict = Dict([region_dem_incoming[i,1],region_dem_incoming[i,2]] => 1 for i in 1:size(region_dem_incoming,1))
param_reg_exc_dict = Dict(reg_exc[i+1,1:2] => reg_exc[i+1,3] for i in 1:size(reg_exc,1)-1)
param_dem_dict = Dict([dem[1,i],j] => dem[j+1,i] for i in 1:size(dem,2), j in 1:size(temporal_resolution,1))
param_cap_dict = Dict(generators[i] => gen_data[i+1,5] for i in 1 : size(generators,1))
param_firm_cap_dict = Dict(generators[i] => 1 for i in 1 : size(generators,1) if gen_data[i+1,10] == "firm")
param_stable_load_dict = Dict(generators[i] => stable_load[i] for i in 1:size(stable_load,1))
param_startup_cost_dict = Dict(generators[i] => startup_cost[i] for i in 1:size(startup_cost,1))
param_gen_cluster_dict = Dict(generators[i] => clusters[i] for i in 1:size(clusters,1))

k = 1
param_capacity_factor_dict = Dict()
for i in 1:size(generators,1)
  if gen_data[i+1,10] != "firm"
   for j in  1:size(temporal_resolution,1)
    param_capacity_factor_dict[generators[i],j] = capacity_factor[k+1,j+2]
  end
 global k += 1
   end
end

param_ramp_up_dict = Dict(generators[i] => ramp[i] for i in 1:size(ramp,1))
param_ramp_down_dict = Dict(generators[i] => ramp[i] for i in 1:size(ramp,1))

power_system = Model(with_optimizer(Gurobi.Optimizer))

@variable(power_system,output[r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution;
         haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i])  && haskey(param_fuel_gen_dict,[f,i])] >=0)

@variable(power_system,com_st[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)] >= 0,Int)
@variable(power_system,start_event[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)] >= 0,Int)
@variable(power_system,shut_event[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)] >= 0,Int)

@constraint(power_system,demand_supply[j in demand_centers,k in temporal_resolution],
            sum(output[r,f,i,j,k] * (haskey(param_firm_cap_dict,i) ? .95 : 1) * (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) for i in generators, r in regions, f in fuels
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])) >= param_dem_dict[[j,k]])

#@constraint(power_system,demand_supply[j in demand_centers,k in temporal_resolution],
            sum(output[r,f,i,j,k] * (haskey(param_reg_dem_dict,[r,j]) ? .95 : .9) for i in generators, r in regions, f in fuels
         if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])) >= param_dem_dict[[j,k]])


@constraint(power_system,capacity_supply[i in generators, k in temporal_resolution],
            sum(output[r,f,i,j,k] for r in regions , j in demand_centers, f in fuels
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                <= (haskey(param_firm_cap_dict,i) ? 0.95 * param_cap_dict[i]*com_st[i,k] : param_cap_dict[i]*param_capacity_factor_dict[i,k]))

#@constraint(power_system,capacity_supply[i in generators, k in temporal_resolution],
            sum(output[r,f,i,j,k] for r in regions , j in demand_centers, f in fuels
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            <= (haskey(param_firm_cap_dict,i) ? 0.85 * param_cap_dict[i] : param_cap_dict[i]*param_capacity_factor_dict[i,k]))

@constraint(power_system,cluster_units[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
            com_st[i,k] <= param_gen_cluster_dict[i])

@constraint(power_system,stable_load[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i)],
            sum(output[r,f,i,j,k] for r in regions , j in demand_centers, f in fuels
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                >= param_stable_load_dict[i] * com_st[i,k] )

@constraint(power_system,start_event_cons[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i) && k != 1],
            start_event[i,k] - shut_event[i,k] == com_st[i,k] - com_st[i,k-1])

@constraint(power_system,start_event_cons_boundary[i in generators, k in temporal_resolution; haskey(param_firm_cap_dict,i) && k == 1],
                        start_event[i,k] == com_st[i,k])

@constraint(power_system,ramp_down[i in generators,k in temporal_resolution; haskey(param_firm_cap_dict,i) && k != temporal_resolution[length(temporal_resolution)]],
            sum(output[r,f,i,j,k] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            -sum(output[r,f,i,j,k+1] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                <= com_st[i,k]*param_ramp_down_dict[i] + param_stable_load_dict[i]*shut_event[i,k])


@constraint(power_system,ramp_up[i in generators,k in temporal_resolution; haskey(param_firm_cap_dict,i) && k != temporal_resolution[length(temporal_resolution)]],
            sum(output[r,f,i,j,k] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
            -sum(output[r,f,i,j,k+1] for r in regions, f in fuels, j in demand_centers
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                >= -com_st[i,k]*param_ramp_up_dict[i] - param_stable_load_dict[i]*start_event[i,k])

@constraint(power_system,transmission_constraint[r in regions, p in regions, k in temporal_resolution; haskey(param_reg_exc_dict,[r,p])],
            sum(output[r,f,i,j,k] for i in generators, f in fuels, j in demand_centers
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_dem_dict,[p,j]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])) <= param_reg_exc_dict[[r,p]])

@objective(power_system, Min, sum(output[r,f,i,j,k]*param_gen_cost_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))
                + sum(start_event[i,k]*param_startup_cost_dict[i] for r in regions, f in fuels, i in generators, k in temporal_resolution
                if haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_firm_cap_dict,i)))

#@objective(power_system, Min, sum(output[r,f,i,j,k]*param_gen_cost_dict[i] for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
            if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])))

optimize!(power_system)

#termination_status(power_system)
#objective_value(power_system)
#value.(output)


#print(power_system)

# writing the model to a file
f = open("power_system.lp", "w")
print(f, power_system)
close(f)

f = open("constraint.lp", "w")
print(f, capacity_supply["I.P. CCPP1.0",1])
close(f)

# writing the results to a file
f = open("results.lp", "w")
print(f, value.(output))
close(f)

# writing results to excel by first converting it to a DataFrame
result_output_df = DataFrame(Region = String[], Fuel = String[], Generator = String[], Demand_center = String[], Time = Int[], Generation = Float64[])
result_start_event_df = DataFrame(Generator = String[], Time = Int[], start_event = Int[])
result_com_st_df = DataFrame(Generator = String[], Time = Int[], Committed = Int[])

for r in regions, f in fuels, i in generators, j in demand_centers, k in temporal_resolution
    if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i])
        push!(result_output_df,(r,f,i,j,k,value(output[r,f,i,j,k])))
        #value(output[i,j,k])
    end
end

result_output_fuel_time_df = DataFrame(Fuel = String[], Time = Int[], Generation = Float64[])
for f in fuels, k in temporal_resolution
    push!(result_output_fuel_time_df,(f,k,sum(value(output[r,f,i,j,k]) for r in regions, i in generators, j in demand_centers
     if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]))))
end


result_output_transmission_df = DataFrame(From = String[], To = String[], Generation = Float64[])
for r in regions, p in regions
    if haskey(param_reg_exc_dict,[r,p])
    push!(result_output_transmission_df,(r,p,sum(value(output[r,f,i,j,k]) for i in generators, j in demand_centers, k in temporal_resolution, f in fuels
     if haskey(param_dem_incoming_dict,[j,r]) && haskey(param_reg_gen_dict,[r,i]) && haskey(param_fuel_gen_dict,[f,i]) && haskey(param_reg_dem_dict,[p,j]))))
end
end

for i in generators, k in temporal_resolution
    if haskey(param_firm_cap_dict,i)
        push!(result_start_event_df,(i,k,value(start_event[i,k])))
        push!(result_com_st_df,(i,k,value(com_st[i,k])))
    end
end

XLSX.writetable("result_output_fuel_time_df.xlsx", collect(DataFrames.eachcol(result_output_fuel_time_df)), DataFrames.names(result_output_fuel_time_df))
XLSX.writetable("result_output_transmission_df.xlsx", collect(DataFrames.eachcol(result_output_transmission_df)), DataFrames.names(result_output_transmission_df))
XLSX.writetable("result_output_df.xlsx", collect(DataFrames.eachcol(result_output_df)), DataFrames.names(result_output_df))
XLSX.writetable("result_start_event_df.xlsx", collect(DataFrames.eachcol(result_start_event_df)), DataFrames.names(result_start_event_df))
XLSX.writetable("result_com_st_df.xlsx", collect(DataFrames.eachcol(result_com_st_df)), DataFrames.names(result_com_st_df))
