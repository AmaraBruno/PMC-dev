#=
    Programa: ScriptPesquisaExpanMercadoCompetitivo
    Grupo de trabalho: PSR - setor de modelos
    
    Descrição:

    Programa criado para implementação de um modelo de otimização para estudo
    da dinâmica do mercado de energia elétrica brasileiro em ambiente totalmente
    competitivo. Neste exemplo, a demanda do modelo é variada e são observados os
    resultados de saída do modelo.
=#


#Importando bibliotecas externas 
using Gurobi, JuMP, PiecewiseLinearOpt, Cbc;

#Importando biblioteca local
include("PMCIO.jl")
include("PMC.jl")

using .PMCIO
using .PMC

####################################################################################################
#Defining parameters
####################################################################################################

#Setting base path to case's paths
toml_file_path="D:\\PMC-dev\\config\\io_config.toml"

#Creating object
io_obj=PMCIO.IOInfo(toml_file_path)
                    
####################################################################################################
#Getting data from input files
####################################################################################################

#Get data
sel_case=10;
num_cases, total_num_plants, num_stages, num_scen, spot_prices, gen_per_plant, pa_per_plant, total_pa_per_case=PMCIO.get_data(io_obj, selected_case=sel_case);

####################################################################################################
#Building PMC model
####################################################################################################

#Model study parameters
demand=11.8;
dem_risk_av_factor=0.7; #lambda_d
gen_risk_av_factor=0.9; #lambda_g
M=10.0e10;
alpha=0.2;
gen_costs=[3.0, 3.0, 3.0, 3.0, 3.0];

#Create PMC object
pmc_model=PMC.PMCModel(demand, 
                       dem_risk_av_factor, 
                       gen_risk_av_factor, 
                       alpha, gen_costs, 
                       gen_per_plant, 
                       pa_per_plant, 
                       total_pa_per_case, 
                       spot_prices, 
                       Cbc.Optimizer,
                       lin_interp=true,
                       print_model=false)

#Show some attributes
println("")
println("PMC demand: $(pmc_model.demand)")
println("PMC number of stages: $(pmc_model.num_stages)")
println("PMC number of scenarios: $(pmc_model.num_scen)")
println("PMC number of plants: $(pmc_model.num_plants)\n")

println("PMC phisical assurances:")
for i=1:length(pmc_model.pa_per_plant)
    println("$(round(pmc_model.pa_per_plant[i], digits=2))    ")
end
println("")

println("PMC spot prices:")
for i=1:size(pmc_model.spot_prices, 2)
    for j=1:size(pmc_model.spot_prices, 3)
        print("$(round(pmc_model.spot_prices[sel_case, i, j], digits=2))    ")
    end
    println("")
end
println("")

####################################################################################################
#Run PMC model and show results
####################################################################################################

demand=[1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 18.0, 21.0]; #lambda_d

#Initialize matrix to store results
VPl_per_plant=zeros(length(demand), total_num_plants);
dispatch_per_plant=zeros(length(demand), total_num_plants);
total_dispatch=zeros(length(demand));
optimal_cost=zeros(length(demand));

#For each demand risk aversion factor, solve the model
i=1
for dem in demand

    #Determine model parameters
    global pmc_model=PMC.PMCModel(dem, 
                                  dem_risk_av_factor, 
                                  gen_risk_av_factor, 
                                  alpha, gen_costs, 
                                  gen_per_plant, 
                                  pa_per_plant, 
                                  total_pa_per_case, 
                                  spot_prices, 
                                  Cbc.Optimizer,
                                  lin_interp=true)

    #Build JuMP model
    GF, Q, Q_plant, hab_plant, suborno, vpl_approx_plant, phi_approx=PMC.build_pmc_model(pmc_model)

    #Run model
    PMC.run_model(pmc_model)

    for j=1:(total_num_plants)
        VPl_per_plant[i,j]=value(vpl_approx_plant[j])
        dispatch_per_plant[i,j]=value(Q_plant[j])
    end

    total_dispatch[i]=value(Q)
    optimal_cost[i]=objective_value(pmc_model.model)

    global i+=1;
end

#Data labels
demand_axis_label="Demand [MW]";
total_dispatch_label=["Hired energy"];
opt_cost_label=["Optimal cost"];
VPl_per_plant_label=["Plant1", "Plant2", "Plant3", "Plant4", "Plant5"];
dispatch_per_plant_label=["Plant1", "Plant2", "Plant3", "Plant4", "Plant5"];

#Creating output
PMCIO.create_output(demand, VPl_per_plant, io_obj, x_label=demand_axis_label, y_label=VPl_per_plant_label, file_name="dem_vpl.csv", output_name="Demanda X Valor presente líquido")
PMCIO.create_output(demand, dispatch_per_plant, io_obj, x_label=demand_axis_label, y_label=dispatch_per_plant_label, file_name="dem_Q_plant.csv", output_name="Demanda X Despacho por planta")
PMCIO.create_output(demand, total_dispatch, io_obj, x_label=demand_axis_label, y_label=total_dispatch_label, file_name="dem_Q.csv", output_name="Demanda X Despacho total")
PMCIO.create_output(demand, optimal_cost, io_obj, x_label=demand_axis_label, y_label=opt_cost_label, file_name="dem_opt_cost.csv", output_name="Demanda X Custo ótimo")
