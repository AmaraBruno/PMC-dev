#=
    Arquivo: PMC

    Grupo de trabalho: PSR - setor de modelos
    
    Descrição:

    Programa criado para implementação de um modelo de otimização para estudo
    da dinâmica do mercado de energia elétrica brasileiro em ambiente totalmente
    competitivo.


=#

module PMC

using Base: Float64
using CSV, DataFrames, JuMP, PiecewiseLinearOpt;

#Exporting functions
export read_data_from_cases
export phi
export vpl

#Defining module-level variables
epsilon=0.5

#Object model
mutable struct PMCModel

    #Study configuration parameters
    demand::Float64
    dem_aversion_risk_factor::Float64 #lambda_d
    gen_aversion_risk_factor::Float64 #lambda_g
    alpha::Float64;
    gen_costs::Vector{Float64}
    lin_interp::Bool 

    #Constant: big M notation
    M::Float64

    #System parameters
    num_stages
    num_scen
    num_plants
    gen_per_plant
    pa_per_plant
    pa_per_case
    spot_prices
    
    #JuMP model variable and optimizer
    optimizer
    model::Model

    #Class constructor
    function PMCModel(demand::Float64,
                      lambda_dem::Float64, 
                      lambda_gen::Float64, 
                      alpha::Float64, 
                      gen_costs::Vector{Float64},
                      gen_per_plant,
                      pa_per_plant,
                      pa_per_case,
                      spot_prices,
                      optimizer;
                      lin_interp::Bool = true)

        #Constant
        M=10.0e10

        #Number of stages and scenarios
        num_stages=size(gen_per_plant, 1);
        num_scen=size(gen_per_plant, 2);

        #Number of plants
        num_plants=size(gen_per_plant, 3)

        #Creating model variable
        model=Model(optimizer)

        #Allocating
        new(demand, 
            lambda_dem, 
            lambda_gen, 
            alpha, 
            gen_costs,
            lin_interp,
            M, 
            num_stages,
            num_scen,
            num_plants,
            gen_per_plant,
            pa_per_plant,
            pa_per_case,
            spot_prices,
            optimizer,
            model)

    end
    
end

function build_pmc_model(pmc_model::PMCModel)

    #Aux
    num_plants=pmc_model.num_plants

    #Adding variables to the model
    @variable(pmc_model.model, hab_plant[1:num_plants], Bin)
    @variable(pmc_model.model, suborno[1:num_plants]>=0.0)
    @variable(pmc_model.model, GF>=0.0)
    @variable(pmc_model.model, Q>=0.0)
    @variable(pmc_model.model, Q_plant[1:num_plants]>=0.0)

    #Adding linear constraints
    @constraint(pmc_model.model, GF==sum(pmc_model.pa_per_plant[i]*hab_plant[i] for i=1:num_plants)) #Soma das GFi viáveis igual a GF
    @constraint(pmc_model.model, Q==sum(Q_plant)) #Soma dos Qi de cada usina igual a Q
    @constraint(pmc_model.model, ger_inviavel[i=1:num_plants], Q_plant[i]<=pmc_model.demand*hab_plant[i]) #Geradores inviáveis não podem ser contratados

    #Defining breakpoints for piecewise approx. for vpl and phi functions
    GF_brkpts=pmc_model.pa_per_case

    #Os montantes contratados não precisam ser restrigindos da mesma forma que a garantia física total
    Q_brkpts=0.0:0.5:35.0;
    Qusina_brkpts=0.0:0.5:35.0;

    #Non-linear constraint for vpl
    if (pmc_model.lin_interp) #If linear interpolation is enabled, UnionJack pattern is not used (Package can evaluate nl function in intermediate points)
        vpl_approx_plant=Array{VariableRef, 1}(undef, num_plants);
        for i=1:num_plants
            vpl_pwl_aux=BivariatePWLFunction(GF_brkpts, Qusina_brkpts, (gf, qi) -> PMC.vpl(gf, qi, pmc_model.gen_costs[i],
            pmc_model.gen_per_plant[:,:,i], pmc_model.pa_per_case, pmc_model.spot_prices, pmc_model.gen_aversion_risk_factor, pmc_model.alpha,
            lin_interp=pmc_model.lin_interp))

            vpl_approx_plant[i]=piecewiselinear(pmc_model.model, GF, Q_plant[i], vpl_pwl_aux)
            @constraint(pmc_model.model, vpl_approx_plant[i]+suborno[i] >= -pmc_model.M*(1.0-hab_plant[i]))
        end

        #Approximating demand expenses
        phi_pwl_aux=BivariatePWLFunction(GF_brkpts, Q_brkpts, (gf, q) -> PMC.phi(gf, q, pmc_model.demand,
        pmc_model.pa_per_case, pmc_model.spot_prices, pmc_model.dem_aversion_risk_factor, pmc_model.alpha, lin_interp=pmc_model.lin_interp))
    else
        vpl_approx_plant=Array{VariableRef, 1}(undef, num_plants);
        for i=1:num_plants
            vpl_pwl_aux=BivariatePWLFunction(GF_brkpts, Qusina_brkpts, (gf, qi) -> PMC.vpl(gf, qi, pmc_model.gen_costs[i],
            pmc_model.gen_per_plant[:,:,i], pmc_model.pa_per_case, pmc_model.spot_prices, pmc_model.gen_aversion_risk_factor, pmc_model.alpha,
            lin_interp=pmc_model.lin_interp), pattern=:UnionJack)

            vpl_approx_plant[i]=piecewiselinear(pmc_model.model, GF, Q_plant[i], vpl_pwl_aux)
            @constraint(pmc_model.model, vpl_approx_plant[i]+suborno[i] >= -pmc_model.M*(1.0-hab_plant[i]))
        end

        #Approximating demand expenses
        phi_pwl_aux=BivariatePWLFunction(GF_brkpts, Q_brkpts, (gf, q) -> PMC.phi(gf, q, pmc_model.demand,
        pmc_model.pa_per_case, pmc_model.spot_prices, pmc_model.dem_aversion_risk_factor, pmc_model.alpha, lin_interp=pmc_model.lin_interp), pattern=:UnionJack)
    end

    #Continuation of demand expenses approx.
    phi_approx=piecewiselinear(pmc_model.model, GF, Q, phi_pwl_aux)

    #Objective function
    @objective(pmc_model.model, Min, phi_approx+sum(suborno[i] for i=1:num_plants))

    #Return JuMP variables
    return GF, Q, Q_plant, hab_plant, suborno, vpl_approx_plant, phi_approx

end

function run_model(pmc_model::PMCModel)

    #Run PMC model
    optimize!(pmc_model.model)

end

function esperanca(vetor_esp_amost,
                   vetor_fdp= fill(1.0/length(vetor_esp_amost), length(vetor_esp_amost)))     
    #Caso vetores informados tenham dimensões diferentes
    if (!(size(vetor_esp_amost)==size(vetor_fdp)))
        error("Input arguments, sample space and prob. density function vectors, differ on size. Verify arguments.")
    end
    
    #Calculando valor esperado
    valor_esperado=0.0
    for i=1:length(vetor_esp_amost)
        valor_esperado+=vetor_fdp[i]*vetor_esp_amost[i]
    end
    
    return valor_esperado
end

#Definindo função CVaR
function cvar(vetor_esp_amost,
              vetor_fdp = fill(1.0/length(vetor_esp_amost), length(vetor_esp_amost));
              alfa,
              pela_esquerda = true)
    #Caso vetores informados tenham dimensões diferentes
    if (!(size(vetor_esp_amost)==size(vetor_fdp)))
    error("Input arguments, sample space and prob. density function vectors, differ on size. Verify arguments.")
    end
    if (abs(alfa)>1.0)
    error("Alfa parameter is defined between 0.0 and 1.0. Verify argument.")
    end
    
    #Inicializando variáveis para o cálculo
    perda_media=0.0;
    tamanho_vetor_esp_amost=length(vetor_esp_amost)
    if (pela_esquerda)
        #Ordenando vetor do espaço amostral em ordem descrescente
        for i=1:tamanho_vetor_esp_amost
            menor=vetor_esp_amost[i];
            indice_menor=i;
            for j=(i+1):tamanho_vetor_esp_amost
                if (vetor_esp_amost[j]<menor)
                    menor=vetor_esp_amost[j];
                    indice_menor=j;
                end
            end
    
            #Caso o índice seja diferente, temos um valor menor e a troca é feita
            aux=0.0;
            if (indice_menor!=i)
                #Troca do vetor do espaço amostral
                aux=vetor_esp_amost[i];
                vetor_esp_amost[i]=vetor_esp_amost[indice_menor];
                vetor_esp_amost[indice_menor]=aux;
                
                #Troca do vetor de densidade de probabilidade
                aux=vetor_fdp[i];
                vetor_fdp[i]=vetor_fdp[indice_menor];
                vetor_fdp[indice_menor]=aux;
            end
        end
    else
        #Ordenando vetor do espaço amostral em ordem descrescente
        for i=1:tamanho_vetor_esp_amost
            maior=vetor_esp_amost[i];
            indice_maior=i;
            for j=(i+1):tamanho_vetor_esp_amost
                if (vetor_esp_amost[j]>maior)
                    maior=vetor_esp_amost[j];
                    indice_maior=j;
                end
            end
    
            #Caso o índice seja diferente, temos um valor menor e a troca é feita
            aux=0.0;
            if (indice_maior!=i)
                #Troca do vetor do espaço amostral
                aux=vetor_esp_amost[i];
                vetor_esp_amost[i]=vetor_esp_amost[indice_maior];
                vetor_esp_amost[indice_maior]=aux;
                
                #Troca do vetor de densidade de probabilidade
                aux=vetor_fdp[i];
                vetor_fdp[i]=vetor_fdp[indice_maior];
                vetor_fdp[indice_maior]=aux;
            end
        end
    end
    
    #Realizando soma ponderada
    perda_media=0.0;
    soma_prob=0.0;
    for i=1:tamanho_vetor_esp_amost
        if ((vetor_fdp[i]+soma_prob) < alfa)
            perda_media+=vetor_fdp[i]*vetor_esp_amost[i];
            soma_prob+=vetor_fdp[i];
        else
            novo_peso=alfa-soma_prob;
            soma_prob+=novo_peso;
            perda_media=(perda_media+novo_peso*vetor_esp_amost[i])/soma_prob;
            break;
        end
    end

    #Retornando
    return perda_media
end

#Definindo função para cálculo do valor presente líquido estocástico
function vpl_estoc_exct(garantia_fisica_total,
                   montante_contratado,
                   custos_geracao,
                   geracoes_por_etapa,
                   vetor_gf_totais, # Vetor de garantias físicas totais associadas a cada matriz
                   matrizes_precos_spot) #Matrizes de preços spot

    #Obtendo número de etapas para iterar sobre
    num_etapas=size(geracoes_por_etapa, 1);

    #Initializing vector to store result
    valor_presente_liquido_estoc=zeros(size(matrizes_precos_spot, 3)) #Tamanho do vetor é o número de cenários
            
    #Para cada garantia física total, faça:
    entrou_if=false;
    for indice_gf=1:length(vetor_gf_totais)
        #println("$garantia_fisica_total - $(vetor_gf_totais[indice_gf])")
        if  (abs(garantia_fisica_total-vetor_gf_totais[indice_gf])<=epsilon)
            #Entrou no if
            entrou_if=true

            #Pré calculando média dos valores da matriz de preços spot
            media_precos_spot=0.0
            for t=1:num_etapas
                media_precos_spot=media_precos_spot+esperanca(matrizes_precos_spot[indice_gf, t,:])
            end
            media_precos_spot=media_precos_spot/Float64(num_etapas);

            for t=1:num_etapas
                valor_presente_liquido_estoc+=(geracoes_por_etapa[t, :] .- montante_contratado).*matrizes_precos_spot[indice_gf, t,:].+
                (media_precos_spot*montante_contratado).-(geracoes_por_etapa[t, :].*custos_geracao);
            end

            #Saia do loop
            break;
        end
    end

    if (!entrou_if)
        error("VPL ESTOC - não entrou no if do loop")
    end

    #Retornando resultado
    return valor_presente_liquido_estoc
end

#Definindo função para cálculo do valor presente líquido estocástico
function vpl_estoc_lin_interp(garantia_fisica_total,
                              montante_contratado,
                              custos_geracao,
                              geracoes_por_etapa,
                              vetor_gf_totais, # Vetor de garantias físicas totais associadas a cada matriz
                              matrizes_precos_spot) #Matrizes de preços spot

    #Obtendo número de etapas para iterar sobre
    num_etapas=size(geracoes_por_etapa, 1);

    #Linear interpolation of spot prices
    precos_spot=zeros(size(matrizes_precos_spot[1, :, :]))
    entrou_if=false
    for i=1:(length(vetor_gf_totais)-1)
        if ((garantia_fisica_total>=vetor_gf_totais[i])&&(garantia_fisica_total<=vetor_gf_totais[i+1]))
            entrou_if=true
            #println("$(garantia_fisica_total)")
            precos_spot=((matrizes_precos_spot[i+1, :, :]-matrizes_precos_spot[i, :, :])./(vetor_gf_totais[i+1]-vetor_gf_totais[i])).*
            (garantia_fisica_total-vetor_gf_totais[i])+matrizes_precos_spot[i, :, :]
            break; #Exit loop
        end
    end

    if (!entrou_if)
        error("Might have exceded bounds. Verify.")
    end

    #Initializing vector to store result
    valor_presente_liquido_estoc=zeros(size(matrizes_precos_spot, 3)) #Tamanho do vetor é o número de cenários
    
    #Pré calculando média dos valores da matriz de preços spot
    media_precos_spot=0.0
    for t=1:num_etapas
        media_precos_spot=media_precos_spot+esperanca(precos_spot[t, :])
    end
    media_precos_spot=media_precos_spot/Float64(num_etapas);

    for t=1:num_etapas
        valor_presente_liquido_estoc+=(geracoes_por_etapa[t, :] .- montante_contratado).*precos_spot[t, :].+
        (media_precos_spot*montante_contratado).-(geracoes_por_etapa[t, :].*custos_geracao);
    end

    #Retornando resultado
    return valor_presente_liquido_estoc
end

#Definindo valor presente líquido
function vpl(garantia_fisica_total,
             montante_contratado,
             custos_geracao,
             geracoes_por_usina,
             vetor_gf_totais,
             matrizes_precos_spot,
             lambda_geracao,
             alfa;
             lin_interp)

    #Verificando argumentos
    if (!(size(geracoes_por_usina, 2)==size(matrizes_precos_spot, 3)))
       error("Scenario dimension mismatch. Verify arguments.") 
    end
    if (abs(lambda_geracao)>1.0)
        error("Risk aversion parameter is defined between 0.0 and 1.0. Verify argument.")
    end
    
    #Variável auxiliar
    aux=0.0;

    if (lin_interp)
        vpl_estoc=vpl_estoc_lin_interp
    else
        vpl_estoc=vpl_estoc_exct
    end
    
    #Calculando vpl estocástico
    valor_presente_liq_estoc=vpl_estoc(garantia_fisica_total,
                                       montante_contratado,
                                       custos_geracao,
                                       geracoes_por_usina,
                                       vetor_gf_totais,
                                       matrizes_precos_spot)

    #Calculando vpl
    aux=lambda_geracao*esperanca(valor_presente_liq_estoc)+(1.0-lambda_geracao)*cvar(valor_presente_liq_estoc, alfa=alfa)
    #Retornando
    return aux
end
    
#Definindo função phi estocástica
function phi_estoc_exct(garantia_fisica_total,
                   montante_contratado,
                   demanda,
                   vetor_gf_totais, # Vetor de garantias físicas totais associadas a cada matriz
                   matrizes_precos_spot) #Matrizes de preços spot

    #Obtendo número de etapas para iterar sobre
    num_etapas=size(matrizes_precos_spot, 2);
    
    #Initializing variable to store result
    gasto_demanda_estoc=zeros(size(matrizes_precos_spot, 3)) #Tamanho do vetor é o número de cenários
    
    #Para cada garantia física total, faça:
    entrou_if=false;
    for indice_gf=1:length(vetor_gf_totais)
        if  (abs(garantia_fisica_total-vetor_gf_totais[indice_gf])<=epsilon)
            #Entrou no if
            entrou_if=true
    
            #Pré calculando média dos valores da matriz de preços spot
            media_precos_spot=0.0
            for t=1:num_etapas
                media_precos_spot=media_precos_spot+esperanca(matrizes_precos_spot[indice_gf, t, :])
            end
            media_precos_spot=media_precos_spot/Float64(num_etapas);
    
            #Calculando vetor phi_estoc
            for t=1:num_etapas
                gasto_demanda_estoc+=(demanda-montante_contratado).*matrizes_precos_spot[indice_gf, t, :].+media_precos_spot*montante_contratado
            end
    
            #Saia do loop
            break;
        end
    end
    
    #Retornando resultado
    return gasto_demanda_estoc
end

function phi_estoc_lin_interp(garantia_fisica_total,
                              montante_contratado,
                              demanda,
                              vetor_gf_totais,      #Vetor de garantias físicas totais associadas a cada matriz
                              matrizes_precos_spot) #Matrizes de preços spot

    #Obtendo número de etapas para iterar sobre
    num_etapas=size(matrizes_precos_spot, 2);

    #Linear interpolation of spot prices
    precos_spot=zeros(size(matrizes_precos_spot[1, :, :]))
    entrou_if=false
    for i=1:(length(vetor_gf_totais)-1)
        if ((garantia_fisica_total>=vetor_gf_totais[i])&&(garantia_fisica_total<=vetor_gf_totais[i+1]))
            entrou_if=true
            precos_spot=((matrizes_precos_spot[i+1, :, :]-matrizes_precos_spot[i, :, :])./(vetor_gf_totais[i+1]-vetor_gf_totais[i])).*
            (garantia_fisica_total-vetor_gf_totais[i])+matrizes_precos_spot[i, :, :]
            break; #Exit loop
        end
    end

    if (!entrou_if)
        error("Might have exceded bounds. Verify.")
    end

    #Using spot prices calculated from the interpolation
    
    #Pré calculando média dos valores da matriz de preços spot
    media_precos_spot=0.0
    for t=1:num_etapas
        media_precos_spot=media_precos_spot+esperanca(precos_spot[t, :])
    end
    media_precos_spot=media_precos_spot/Float64(num_etapas);

    #Initializing variable to store result
    gasto_demanda_estoc=zeros(size(matrizes_precos_spot, 3)) #Tamanho do vetor é o número de cenários
    
    #Calculando vetor phi_estoc
    for t=1:num_etapas
        gasto_demanda_estoc+=(demanda-montante_contratado).*precos_spot[t, :].+media_precos_spot*montante_contratado
    end

    #Retornando resultado
    return gasto_demanda_estoc
end


#Definindo função phi
function phi(garantia_fisica_total,
             montante_contratado,
             demanda,
             vetor_gf_totais, # Vetor de garantias físicas totais associadas a cada matriz
             matrizes_precos_spot, #Matrizes de preços spot
             lambda_demanda,
             alfa;
             lin_interp)

    #Verificando argumentos
    if (abs(alfa)>1.0)
        error("Alfa parameter is defined between 0.0 and 1.0. Verify argument.")
    end
    if (abs(lambda_demanda)>1.0)
        error("Risk aversion parameter is defined between 0.0 and 1.0. Verify argument.")
    end
    
    if (lin_interp)
        phi_estoc=phi_estoc_lin_interp
    else
        phi_estoc=phi_estoc_exct
    end

    #Variável auxiliar
    aux=0.0;
    
    #Calculando phi estocástico
    gasto_demanda_estoc=phi_estoc(garantia_fisica_total, montante_contratado, demanda, vetor_gf_totais,  matrizes_precos_spot)
    
    #Calculando vpl. Sentido do CVaR é para cima
    aux=lambda_demanda*esperanca(gasto_demanda_estoc)+(1.0-lambda_demanda)*cvar(gasto_demanda_estoc, alfa=alfa, pela_esquerda=false)
    
    #Retornando
    return aux
end

end #Fim módulo