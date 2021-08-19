#Create module
module PMCIO

using DataFrames: Matrix
using CSV, DataFrames, Plots, TOML;

export IOInfo, getData

#Define struct to store info about files
mutable struct IOInfo
    #Paths
    #path_to_cases_dir_file::String
    #path_to_data_files_names::String
    io_config_file_path::String
    output_directory_path::String
    
    #File IO config
    #open_mode::String
    #datasep::String
    initial_line::Int
    index_dict::Dict

    #Struct constructor
    function IOInfo(#p2p_file::String,
                    #path2_data_files::String,
                    toml_file_path)
        #Configuring IO specs
        #open_mode = "r";
        #datasep = ",";

        #Parse toml config file
        content=TOML.parsefile(toml_file_path);
        output_dir=content["pmc_model_output_dir"]["output_path_dir"];
        
        initial_line = 4; #header line position
        index_dict = Dict("cmgdem.csv" => 1, "gerter.csv" => 2, "gerhid.csv" => 3);

        #Construct IOInfo object
        new(toml_file_path,
            output_dir, 
            initial_line,
            index_dict);
            
    end
end

function get_data(stct::IOInfo; selected_case::Int, verbose::Bool=true)

    #IO params
    initial_line=stct.initial_line
    index_dict=stct.index_dict

    #Parse file
    content=TOML.parsefile(stct.io_config_file_path);

    #Get file's content
    files=content["sddp_files"]["file_name"];
    cases=content["sddp_cases"]["paths"];
    
    num_cases=length(cases);
    num_files=length(files);

    #Getting data in files
    data_stct=Matrix{DataFrame}(undef, (num_cases, num_files))
    for i=1:num_cases
    println("Opening case $i of $num_cases - ",round(100.0*i/num_cases), "%:")
    for j=1:num_files
            #Showing file
            file_name=files[j];
            full_path=string(cases[i], "\\", file_name);
            println(full_path);
            
            #Get data
            data_stct[i,index_dict[file_name]]=CSV.read(full_path, DataFrame, header=initial_line); 
    end
    println("")
    end

    #Assuming that number of plants does not change between files and cases
    num_stages=maximum(data_stct[1, 1][:,1])
    num_scen=maximum(data_stct[1, 1][:, 2])
    num_therm_plts=ncol(data_stct[1, index_dict["gerter.csv"]])-3
    num_hyd_plts=ncol(data_stct[1, index_dict["gerhid.csv"]])-3
    total_num_plts=num_therm_plts+num_hyd_plts;

    #Alocating memory
    spot_prices=zeros(num_cases, num_stages, num_scen);
    therm_gen=zeros(num_cases, num_therm_plts, num_stages, num_scen);
    hyd_gen=zeros(num_cases, num_hyd_plts, num_stages, num_scen);
    gen_per_plant=zeros(num_stages, num_scen, total_num_plts);

    #Opening data structures
    for case=1:num_cases
        #Spot prices
        for row in eachrow(data_stct[case, index_dict["cmgdem.csv"]])
            spot_prices[case, row[1], row[2]]=row[4];
        end
        
        #Thermal generation
        for row in eachrow(data_stct[case, index_dict["gerter.csv"]])
            for i=1:num_therm_plts
                therm_gen[case, i, row[1], row[2]]=row[i+3];
            end
        end
        
        #Hydro generation
        for row in eachrow(data_stct[case, index_dict["gerhid.csv"]])
            for i=1:num_hyd_plts
                hyd_gen[case, i, row[1], row[2]]=row[i+3];
            end
        end
    end

    #Calculating generation matrix for each plant
    for i=1:num_therm_plts
        gen_per_plant[:, :, i] =  therm_gen[selected_case, i, :, :];
    end
    for i=(num_therm_plts+1):total_num_plts
        gen_per_plant[:, :, i] =  hyd_gen[selected_case, i-num_therm_plts, :, :];
    end

    #Create variables to store results
    therm_pa = Matrix{Float64}(undef, (num_cases, num_therm_plts))
    hyd_pa = Matrix{Float64}(undef, (num_cases, num_hyd_plts))
    total_pa_per_case = Vector{Float64}(undef, num_cases)
    pa_per_plant=Vector{Float64}(undef, total_num_plts);

    for i=1:num_cases
        #Estimating phisical assurance of the thermal plants
        for j=1:num_therm_plts
            therm_pa[i, j]=sum(therm_gen[i, j, :, :])/(num_stages*num_scen)
        end
        
        #Estimating phisical assurance of the hydroelectric plants
        for j=1:num_hyd_plts
            hyd_pa[i, j]=sum(hyd_gen[i, j, :, :])/(num_stages*num_scen)
        end
    end

    #Concatenating phisical assurances of power plants
    for i=1:(num_therm_plts)
        pa_per_plant[i]=therm_pa[i];
    end

    for i=(num_therm_plts+1):(total_num_plts)
        pa_per_plant[i]=hyd_pa[i-num_therm_plts];
    end

    #Total phisical assurance per case
    for i=1:num_cases
        total_pa_per_case[i]=sum(hyd_pa[i,:])+sum(therm_pa[i, :])
    end

    if (verbose)
        println("Number of cases: $num_cases");
        println("Number of files to read: $num_files");
        println("Number of stages: $num_stages")
        println("Number of scenarios: $num_scen")
        println("Number of thermal plants: $num_therm_plts")
        println("Number of hydro plants: $num_hyd_plts")
        println("Total number of plants: $total_num_plts")
    end

    return num_cases, total_num_plts, num_stages, num_scen, spot_prices, gen_per_plant, pa_per_plant, total_pa_per_case

end
#=
#Functions to operate over the class
function get_data(stct::IOInfo)

    #Get configs
    cases_dirs_path_file=stct.path_to_cases_dir_file
    path_to_data_files_names=stct.path_to_data_files_names

    open_mode=stct.open_mode
    datasep=stct.datasep
    initial_line=stct.initial_line
    index_dict=stct.index_dict

    #Getting cases path
    num_cases=0
    cases_paths=String[] #Initialize vector
    println("Opening file containing paths to the cases:\n$(cases_dirs_path_file)\n")
    cases_path_file=open(cases_dirs_path_file, open_mode)
    for line in readlines(cases_path_file)
    num_cases+=1
    resize!(cases_paths, num_cases)
    cases_paths[num_cases]=line
    end

    #Verifying data
    cases_size=length(cases_paths)
    if (cases_size==0)
        error("No path informed to cases. Verify file.")
    end

    #Showing cases paths
    println("Number of cases: $cases_size")
    println("Cases directories:")
    for case in cases_paths
    println(case) 
    end

    #Reading file containing cases data files to be read
    num_files=0
    file_names_aux=String[]
    data_files_name_file=open(path_to_data_files_names, open_mode)
    for line in readlines(data_files_name_file)
        num_files+=1
        resize!(file_names_aux, num_files*2)
        keyword, file_name = lstrip.(split(line, datasep))
        file_names_aux[2*(num_files-1)+1]=keyword
        file_names_aux[2*(num_files-1)+2]=file_name
    end

    #Initialing matrix to store data
    file_names=["Something" for i=1:num_files, j=1:2]
    for i=1:num_files
        file_names[i, 1]=file_names_aux[2*(i-1)+1]
        file_names[i, 2]=file_names_aux[2*(i-1)+2]
    end
    println("")

    #Printing info got from the file
    println("Number of files to read: $num_files")
    println("Data from $(path_to_data_files_names):")
    for i=1:size(file_names, 1)
        for j=1:size(file_names, 2)
            print(file_names[i, j], " ")
        end
        println("")
    end
    println("")

    #Getting data in files
    data_stct=Matrix{DataFrame}(undef, (num_cases, num_files))
    for i=1:num_cases
    println("Opening files from case $i of $num_cases - ",round(100.0*i/num_cases), "%:")
    for j=1:num_files
            #Showing file
            full_path=string(cases_paths[i], "\\", file_names[j, 2])
            println(full_path)
            
            #Get type of file
            file_type=file_names[j, 1]
            
            #Get data
            data_stct[i,index_dict[file_type]]=CSV.read(full_path, DataFrame, header=initial_line)        
    end
    println("")
    end

    #Assuming that number of plants does not change between files and cases
    num_stages=maximum(data_stct[1, 1][:,1])
    num_scen=maximum(data_stct[1, 1][:, 2])
    num_therm_plts=ncol(data_stct[1, index_dict["TER"]])-3
    num_hyd_plts=ncol(data_stct[1, index_dict["HID"]])-3

    println("Number of stages: $num_stages")
    println("Number of scenarios: $num_scen")
    println("Number of thermal plants: $num_therm_plts")
    println("Number of hydro plants: $num_hyd_plts")

    #Alocating memory
    spot_prices=zeros(num_cases, num_stages, num_scen);
    therm_gen=zeros(num_cases, num_therm_plts, num_stages, num_scen);
    hyd_gen=zeros(num_cases, num_hyd_plts, num_stages, num_scen);

    #Opening data structures
    for case=1:num_cases
        #Spot prices
        for row in eachrow(data_stct[case, index_dict["CMG"]])
            spot_prices[case, row[1], row[2]]=row[4];
        end
        
        #Thermal generation
        for row in eachrow(data_stct[case, index_dict["TER"]])
            for i=1:num_therm_plts
                therm_gen[case, i, row[1], row[2]]=row[i+3];
            end
        end
        
        #Hydro generation
        for row in eachrow(data_stct[case, index_dict["HID"]])
            for i=1:num_hyd_plts
                hyd_gen[case, i, row[1], row[2]]=row[i+3];
            end
        end
    end

    return spot_prices,
           therm_gen,
           hyd_gen,
           num_cases,
           num_stages,
           num_scen, 
           num_therm_plts,
           num_hyd_plts

end

function process_data(therm_gen, hyd_gen)

    #Getting parameters
    num_cases=size(therm_gen, 1)
    num_therm_gen=size(therm_gen, 2)
    num_hyd_gen=size(hyd_gen, 2)
    num_stages=size(therm_gen, 3)
    num_scen=size(therm_gen, 4)

    #Create variables to store results
    therm_pa = Matrix{Float64}(undef, (num_cases, num_therm_gen))
    hyd_pa = Matrix{Float64}(undef, (num_cases, num_hyd_gen))
    total_pa_per_case = Vector{Float64}(undef, num_cases)

    for i=1:num_cases
        #Estimating phisical assurance of the thermal plants
        for j=1:num_therm_gen
            therm_pa[i, j]=sum(therm_gen[i, j, :, :])/(num_stages*num_scen)
        end
        
        #Estimating phisical assurance of the hydroelectric plants
        for j=1:num_hyd_gen
            hyd_pa[i, j]=sum(hyd_gen[i, j, :, :])/(num_stages*num_scen)
        end
    end
    
    #Total phisical assurance per case
    for i=1:num_cases
        total_pa_per_case[i]=sum(hyd_pa[i,:])+sum(therm_pa[i, :])
    end

    return therm_pa, hyd_pa, total_pa_per_case

end

function calc_per_plant(therm_gen, hyd_gen, therm_pa, hyd_pa)

    num_stages=size(therm_gen, 2);
    num_scen=size(therm_gen, 3);
    num_therm=size(therm_gen, 1);
    num_hyd=size(hyd_gen, 1);
    num_plants=num_therm+num_hyd;

    pa_per_plant=Vector{Float64}(undef, num_therm+num_hyd);
    gen_per_plant=zeros(num_stages, num_scen, num_plants);

    for i=1:(num_therm)
        pa_per_plant[i]=therm_pa[i];
    end
    for i=(num_therm+1):(num_therm+num_hyd)
        pa_per_plant[i]=hyd_pa[i-num_therm];
    end
    
    #Calculating generation matrix for each plant
    for i=1:num_therm
        gen_per_plant[:, :, i] =  therm_gen[i, :, :];
    end
    for i=(num_therm+1):num_plants
        gen_per_plant[:, :, i] =  hyd_gen[i-num_therm, :, :];
    end
    
    return pa_per_plant, gen_per_plant;

end
=#
function create_output(x, y, stct::IOInfo;x_label, y_label, file_name, output_name)
    
    #Configure plot
    font_name="times"
    Plots.default(titlefont = (font_name),
                  xguidefont=(font_name),
                  yguidefont=(font_name),
                  legendfont = (font_name), 
                  xtickfont=(font_name), 
                  ytickfont=(font_name))
                  
    #Get columns of matrix y
    num_cols=size(y, 2);
    charts_path=string(stct.output_directory_path, "\\Charts");
    
    if (length(y[1, :])!=num_cols)
        error("Number of y labels and columns do not match")
    end
    
    #If file extension was included, remove it
    if occursin(".csv", file_name)
        file_name=replace(file_name, ".csv"=>"")
    end

    #Initialize DataFrames obj
    df=DataFrame();
    
    #Add column for x values
    pos=1;
    insertcols!(df, pos, x_label => vec(x));
    
    #Add column for y values
    for i=1:num_cols
       pos+=1
       insertcols!(df, pos, y_label[i] => vec(y[:, i]))
    end

    #Verify folder existence
    if (!isdir(stct.output_directory_path))
        mkdir(stct.output_directory_path)
        mkdir(charts_path)
    end

    #If Charts folder does not exist
    if (!isdir(charts_path))
        mkdir(charts_path)
    end
        
    #Save dataframe to CSV file
    CSV.write(string(stct.output_directory_path, "\\", file_name, ".csv"), df)

    #Create chart
    plot(x, 
         y, 
         title=output_name, 
         xlabel=x_label,
         label=reshape(y_label, 1, length(y_label)),
         legend=:outertopright, 
         lw=3)

    savefig(string(charts_path,"\\", output_name, ".png"))
    
end

end #End module IO