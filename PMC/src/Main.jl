#=
    Programa: ScriptPesquisaExpanMercadoCompetitivo
    Grupo de trabalho: PSR - setor de modelos
    
    Descrição:

    Programa criado para implementação de um modelo de otimização para estudo
    da dinâmica do mercado de energia elétrica brasileiro em ambiente totalmente
    competitivo.
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
path_to_path_file="D:\\Dropbox (PSR)\\MeusProgramas\\Julia\\PesquisaEquilibrioMercado\\CodigosJulia\\src\\config_files\\cases_paths.txt"
path_to_data_files="D:\\Dropbox (PSR)\\MeusProgramas\\Julia\\PesquisaEquilibrioMercado\\CodigosJulia\\src\\config_files\\data_files_names.txt"
output_dir="D:\\Dropbox (PSR)\\MeusProgramas\\Julia\\PesquisaEquilibrioMercado\\CodigosJulia\\src\\Results"

#Creating object
io_obj=PMCIO.IOInfo(path_to_path_file,
                    path_to_data_files,
                    output_dir)
                    
####################################################################################################
#Getting data from input files
####################################################################################################

#Get data
spot_prices,therm_gen,hyd_gen,num_cases,num_stages,num_scen,num_therm_plts,num_hyd_plts=PMCIO.get_data(io_obj);

####################################################################################################
#Data processing
####################################################################################################

#Process data
therm_pa, hyd_pa, total_pa_per_case=PMCIO.process_data(therm_gen, hyd_gen)

#Showing Results
println("")
println("Thermal phisical assurances:")
for i=1:size(therm_pa, 1)
    for j=1:size(therm_pa, 2)
        print("$(round(therm_pa[i, j], digits=2))    ")
    end
    println("")
end
println("")

println("Hydro phisical assurances:")
for i=1:size(hyd_pa, 1)
    for j=1:size(hyd_pa, 2)
        print("$(round(hyd_pa[i, j], digits=2))    ")
    end
    println("")
end
println("")

println("Total phisical assurances per case:")
for i=1:length(total_pa_per_case)
    println("Case $i: $(total_pa_per_case[i])")
end

#Revert order to simulate PA and spot prices relationship
spot_prices=reverse(spot_prices, dims=1)

#Calculate pa per plant and generation matrix per plant
sel_case=5
pa_per_plant, gen_per_plant=PMCIO.calc_per_plant(therm_gen[sel_case, :, :, :], hyd_gen[sel_case, :, :, :], therm_pa[sel_case, :], hyd_pa[sel_case, :]);

println("")
println("Phisical assurance per plant")
for i=1:length(pa_per_plant)
    println("$(pa_per_plant[i])")
end
println("")

println("Generation of plant 1")
for i=1:size(gen_per_plant, 1)
    for j=1:size(gen_per_plant, 2)
        print("$(gen_per_plant[i, j, 1])   ")
    end
    println("")
end

####################################################################################################
#Building PMC model
####################################################################################################

#Model study parameters
demand=11.8;
dem_risk_av_factor=0.7; #lambda_d
gen_risk_av_factor=0.9; #lambda_g
M=10.0e10;
alpha=0.2;
gen_costs=[0.8, 0.8, 1.2, 0.0, 0.0];

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
                       lin_interp=true)

#Show some attributes
println("")
println("PMC demand: $(pmc_model.demand)")
println("PMC number of stages: $(pmc_model.num_stages)")
println("PMC number of scenarios: $(pmc_model.num_scen)")
println("PMC number of plants: $(pmc_model.num_plants)\n")

println("PMC spot prices:")
for i=1:size(pmc_model.spot_prices, 2)
    for j=1:size(pmc_model.spot_prices, 3)
        print("$(round(pmc_model.spot_prices[5, i, j], digits=2))    ")
    end
    println("")
end
println("")

####################################################################################################
#Run PMC model and show results
####################################################################################################

demand=[1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 18.0, 21.0]; #lambda_d

#Initialize matrix to store results
VPl_per_plant=zeros(length(demand), num_hyd_plts+num_therm_plts);
dispatch_per_plant=zeros(length(demand), num_hyd_plts+num_therm_plts);
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

    for j=1:(num_hyd_plts+num_therm_plts)
        VPl_per_plant[i,j]=value(vpl_approx_plant[j])
        dispatch_per_plant[i,j]=value(Q_plant[j])
    end

    total_dispatch[i]=value(Q)
    optimal_cost[i]=objective_value(pmc_model.model)

    global i+=1;
end

#Data labels
av_dem_factor_label="Demand aversion";
total_dispatch_label=["Hired energy"];
opt_cost_label=["Optimal cost"];
VPl_per_plant_label=["Plant1", "Plant2", "Plant3", "Plant4", "Plant5"];
dispatch_per_plant_label=["Plant1", "Plant2", "Plant3", "Plant4", "Plant5"];

#Creating output
PMCIO.create_output(demand, VPl_per_plant, io_obj, x_label=av_dem_factor_label, y_label=VPl_per_plant_label, file_name="dem_vpl.csv", output_name="Demanda X Valor presente líquido")
PMCIO.create_output(demand, dispatch_per_plant, io_obj, x_label=av_dem_factor_label, y_label=dispatch_per_plant_label, file_name="dem_Q_plant.csv", output_name="Demanda X Despacho por planta")
PMCIO.create_output(demand, total_dispatch, io_obj, x_label=av_dem_factor_label, y_label=total_dispatch_label, file_name="dem_Q.csv", output_name="Demanda X Despacho total")
PMCIO.create_output(demand, optimal_cost, io_obj, x_label=av_dem_factor_label, y_label=opt_cost_label, file_name="dem_opt_cost.csv", output_name="Demanda X Custo ótimo")
