function [ModelResults, solveroutput] = DAC_Standalone_6(Solver_Options, Settings, InputData, RH_Vector, T_Vector, Input_DAC, SaveDetails, upperlimit)
    tic
    %------ SET DATA
    % Efficiency of ohmig resistance (transformation from electrical to
    % thermal energy
    eta_elth        = InputData.eta_elth;
    
    % Trade off between electricity and heat?
    TradeOff = Settings.ElThTradeOff;
    
    % Price Data
    P_el            = InputData.EconomicData.CarrierPrices.Electricity;
    P_th            = InputData.EconomicData.CarrierPrices.Heat;
    % Interest rate
    r               = InputData.EconomicData.DAC.InterestRate;
    % Lifetime of the investment
    LT              = InputData.EconomicData.DAC.Lifetime;
    % CO2 demand
    Demand          = InputData.Demand;
    % CO2 Factors
    CO2FactorEl     = InputData.EmissionFactors.Electricity;
    CO2FactorTh     = InputData.EmissionFactors.Heat;
    % Size constraints
    Model.size_min  = Settings.Size_min;
    Model.size_max  = Settings.Size_max;
    % Constant operating parameters
    ConstantOP = Settings.ConstantOP;
    RHadapt = Settings.RHadapt;
    RHdiscret = Settings.RHdiscret;
    
    
    % Time settings
    N               = Settings.N;  % number of days per year 
    D               = Settings.D;  % number of design hours
    K               = Settings.K;  % number of daily intervals
    
    if ConstantOP
        nrOS  = unique(Input_DAC.OSIdentifier1);
        nrOS  = size(nrOS,1);
    else
        % Number of breakpoints on the piecewise approximated performance
        % curves
        nr_bp           = Settings.nr_bp;
        % Number of operational strategies
        nrOS  = unique(Input_DAC.OS);
        nrOS  = size(nrOS,1);
        % One or multiple OS?
        OS = Settings.OS;
        if RHadapt
            OS = 2;
        end
        nrPieces = nr_bp-1;
    end
    
    
    
    %------ K-MEANS CLUSTERING
    % Reshaping Vectors
    RH              = reshape(RH_Vector, [K, N])';
    T               = reshape(T_Vector, [K, N])';
    P_el            = reshape(P_el, [K, N])';
    P_th            = reshape(P_th, [K, N])';
    Demand          = reshape(Demand, [K, N])';
    CO2FactorEl     = reshape(CO2FactorEl, [K, N])';
    CO2FactorTh     = reshape(CO2FactorTh, [K, N])';

    % Vector to perform k-means clustering on
    x = [RH T P_el P_th Demand CO2FactorEl CO2FactorTh];

    % Create a starting vector
    start = floor(linspace(1, N, D));
    % Perform clustering
    [Model.ClusteredData.idx, c] = kmeans(x,D,'start', x(start,:));

    % Retrieve data from k-means matrix
    i=0;
    Model.ClusteredData.RH          = c(:,i*K+1:(i+1)*K);
    i=i+1;
    Model.ClusteredData.T           = c(:,i*K+1:(i+1)*K);
    i=i+1;
    Model.ClusteredData.P_el        = c(:,i*K+1:(i+1)*K);
    i=i+1;
    Model.ClusteredData.P_th        = c(:,i*K+1:(i+1)*K);
    i=i+1;
    Model.ClusteredData.Demand      = c(:,i*K+1:(i+1)*K);
    i=i+1;
    Model.ClusteredData.CO2FactorEl = c(:,i*K+1:(i+1)*K);
    i=i+1;
    Model.ClusteredData.CO2FactorTh = c(:,i*K+1:(i+1)*K);

    %% ------ PREPROCESSING
    % unit adjustments
    Input_DAC.E_tot   = Input_DAC.CO2_Out .* Input_DAC.E_tot /3.6; % in kwh/h
    Input_DAC.E_el    = Input_DAC.CO2_Out .* Input_DAC.E_el /3.6; % in kwh/h
    Input_DAC.Water_Out = Input_DAC.Water_Out .* Input_DAC.CO2_Out/1000; % in t/h
    Input_DAC.CO2_Out = Input_DAC.CO2_Out/1000; % in t/h
    if ConstantOP
        if RHadapt
            load('AirinletCooling.mat')
            temp_ind = (Model.ClusteredData.T >= 15 & Model.ClusteredData.RH <= 60);
            RH_adapt = Model.ClusteredData.RH;
            RH_adapt(temp_ind) = 60;
            F = scatteredInterpolant(AirinletCoolingS1.T_amb, AirinletCoolingS1.RH_amb, AirinletCoolingS1.RH_DAC, AirinletCoolingS1.T_DAC, 'linear', 'none');
            T_adapt = F(Model.ClusteredData.T, Model.ClusteredData.RH, RH_adapt);
            F = scatteredInterpolant(AirinletCoolingS1.T_amb, AirinletCoolingS1.RH_amb, AirinletCoolingS1.RH_DAC, AirinletCoolingS1.M_water_kgs, 'linear', 'none');
            Water_cons(:, :, i) = F(Model.ClusteredData.T, Model.ClusteredData.RH, RH_adapt)/10;
            Model.ClusteredData.T(temp_ind) = T_adapt(temp_ind);
            Model.ClusteredData.RH(temp_ind) = RH_adapt(temp_ind);
            [alpha, epsilon_el, epsilon_th, ~] = DACOutput_4_ConsPar(Model.ClusteredData.T, Model.ClusteredData.RH, Input_DAC);
        else
            [alpha, epsilon_el, epsilon_th, ~] = DACOutput_4_ConsPar(Model.ClusteredData.T, Model.ClusteredData.RH, Input_DAC);
        end
        
    else
            if RHadapt
                nrOS = RHdiscret;
            else
                nrOS = 1;
            end
            b_i_1 = zeros(D, K, nrOS, nrPieces);
            b_i   = zeros(D, K, nrOS, nrPieces);
            alpha = zeros(D, K, nrOS, nrPieces);
            beta  = zeros(D, K, nrOS, nrPieces);
            a_i_1 = zeros(D, K, nrOS, nrPieces);
            a_i   = zeros(D, K, nrOS, nrPieces);
            gamma = zeros(D, K, nrOS, nrPieces);
            delta = zeros(D, K, nrOS, nrPieces);
            temp_data = Input_DAC(Input_DAC.OS == OS,:);
            
            if RHadapt
                load('AirinletCooling.mat')
                RH_choices = zeros(D, K, nrOS);
                T_choices = zeros(D, K, nrOS);
                Water_cons = zeros(D, K, nrOS);
                RH_disc = linspace(0,100, nrOS);
                [b_i_1RH(:, :, 1, :), b_iRH(:, :, 1, :), alphaRH(:, :, 1, :), betaRH(:, :, 1, :), a_i_1RH(:, :, 1, :), a_iRH(:, :, 1, :), gammaRH(:, :, 1, :), deltaRH(:, :, 1, :), ~, ~, ~, ~] =...
                        DACOutput_4(Model.ClusteredData.T, Model.ClusteredData.RH, nr_bp, temp_data);
                for i = 1:nrOS
                    RH_adapt1 = ones(D, K) .* RH_disc(i);
                    if i == nrOS
                        RH_adapt2 = RH_adapt1+1;
                    else
                        RH_adapt2 = ones(D, K) .* RH_disc(i+1);
                    end
                    check1 = RH_adapt1 < Model.ClusteredData.RH;
                    check2 = (RH_adapt1 <= Model.ClusteredData.RH) & (Model.ClusteredData.RH < RH_adapt2);
                    %We do not assume that humidity can be increased for
                    %temperatures below 15°C, productivity will be set to
                    %zero then
                    F = scatteredInterpolant(AirinletCoolingS1.T_amb, AirinletCoolingS1.RH_amb, AirinletCoolingS1.RH_DAC, AirinletCoolingS1.T_DAC, 'linear', 'none');
                    T_adapt = F(Model.ClusteredData.T, Model.ClusteredData.RH, RH_adapt1);
                    F = scatteredInterpolant(AirinletCoolingS1.T_amb, AirinletCoolingS1.RH_amb, AirinletCoolingS1.RH_DAC, AirinletCoolingS1.M_water_kgs, 'linear', 'none');
                    Water_cons(:, :, i) = F(Model.ClusteredData.T, Model.ClusteredData.RH, RH_adapt1)/10;
                    check3 = isnan(T_adapt);
                    [b_i_1(:, :, i, :), b_i(:, :, i, :), alpha(:, :, i, :), beta(:, :, i, :), a_i_1(:, :, i, :), a_i(:, :, i, :), gamma(:, :, i, :), delta(:, :, i, :), ~, ~, ~, ~] =...
                        DACOutput_4(T_adapt, RH_adapt1, nr_bp, temp_data);
                    temp_RH = RH_adapt1;
                    temp_RH(check2) = Model.ClusteredData.RH(check2);
                    RH_choices(:,:,i) = temp_RH;
                    
                    
                    temp_T = T_adapt;
                    temp_T(check2) = Model.ClusteredData.T(check2);
                    T_choices(:,:,i) = temp_T;
                    for j=1:nrPieces
                        temp_b_i_1 = b_i_1(:, :, i, j);
                        temp_b_i_1RH = b_i_1RH(:, :, 1, j);
                        % Set performance to zero if RHadapt is below RHactual
                        temp_b_i_1(check1)=0;
                        % Set performance to zero if Tadapt yields extrapolation
                        temp_b_i_1(check3)=0;
                        % Set performance to performance of DAC without RH changes 
                        temp_b_i_1(check2)=temp_b_i_1RH(check2);
                        b_i_1(:, :, i, j) = temp_b_i_1;
                        
                        temp_b_i = b_i(:, :, i, j);
                        temp_b_iRH = b_iRH(:, :, 1, j);
                        temp_b_i(check1)=0;
                        temp_b_i(check3)=0;
                        temp_b_i(check2)=temp_b_iRH(check2);
                        b_i(:, :, i, j) = temp_b_i;
                        
                        temp_alpha = alpha(:, :, i, j);
                        temp_alphaRH = alphaRH(:, :, 1, j);
                        temp_alpha(check1)=0;
                        temp_alpha(check3)=0;
                        temp_alpha(check2)=temp_alphaRH(check2);
                        alpha(:, :, i, j) = temp_alpha;
                        
                        temp_beta = beta(:, :, i, j);
                        temp_betaRH = betaRH(:, :, 1, j);
                        temp_beta(check1)=0;
                        temp_beta(check3)=0;
                        temp_beta(check2)=temp_betaRH(check2);
                        beta(:, :, i, j) = temp_beta;
                        
                        temp_a_i_1 = a_i_1(:, :, i, j);
                        temp_a_i_1RH = a_i_1RH(:, :, 1, j);
                        temp_a_i_1(check1)=0;
                        temp_a_i_1(check3)=0;
                        temp_a_i_1(check2)=temp_a_i_1RH(check2);
                        a_i_1(:, :, i, j) = temp_a_i_1;
                        
                        temp_a_i = a_i(:, :, i, j);
                        temp_a_iRH = a_iRH(:, :, 1, j);
                        temp_a_i(check1)=0;
                        temp_a_i(check3)=0;
                        temp_a_i(check2)=temp_a_iRH(check2);
                        a_i(:, :, i, j) = temp_a_i;
                        
                        temp_gamma = gamma(:, :, i, j);
                        temp_gammaRH = gammaRH(:, :, 1, j);
                        temp_gamma(check1)=0;
                        temp_gamma(check3)=0;
                        temp_gamma(check2)=temp_gammaRH(check2);
                        gamma(:, :, i, j) = temp_gamma;
                        
                        temp_delta = delta(:, :, i, j);
                        temp_deltaRH = deltaRH(:, :, 1, j);
                        temp_delta(check1)=0;
                        temp_delta(check3)=0;
                        temp_delta(check2)=temp_deltaRH(check2);
                        delta(:, :, i, j) = temp_delta;
                        
                    end
                end
                if nrPieces == 1
                    ind_temp = squeeze(sum(sum(alpha,1),2)>0);
                    ind_temp2 = ind_temp .* [1:size(ind_temp,1)]';
                    minOS = min(ind_temp2(ind_temp2>0));
                    maxOS = max(ind_temp2);

                    b_i_1   = b_i_1(:, :, minOS:maxOS);
                    b_i     = b_i(:, :, minOS:maxOS);
                    alpha   = alpha(:, :, minOS:maxOS);
                    beta    = beta(:, :, minOS:maxOS);
                    a_i_1   = a_i_1(:, :, minOS:maxOS);
                    a_i     = a_i(:, :, minOS:maxOS);
                    gamma   = gamma(:, :, minOS:maxOS);
                    delta   = delta(:, :, minOS:maxOS);
                    RH_choices = RH_choices(:, :, minOS:maxOS);
                    T_choices = T_choices(:, :, minOS:maxOS);
                    Water_cons = Water_cons(:, :, minOS:maxOS);
                    nrOS = sum(ind_temp);
                end
            else
                [b_i_1(:, :, 1, :), b_i(:, :, 1, :), alpha(:, :, 1, :), beta(:, :, 1, :), a_i_1(:, :, 1, :), a_i(:, :, 1, :), gamma(:, :, 1, :), delta(:, :, 1, :), ~, ~, ~, ~] =...
                    DACOutput_4(Model.ClusteredData.T, Model.ClusteredData.RH, nr_bp, temp_data);
                Model.OS = OS;
                
                
                below5 = 0;
                if below5 == 1
                    b_i_1(Model.ClusteredData.T < 5)    = 0;
                    b_i(Model.ClusteredData.T < 5)      = 0;
                    alpha(Model.ClusteredData.T < 5)    = 0;
                    beta(Model.ClusteredData.T < 5)     = 0;
                    a_i_1(Model.ClusteredData.T < 5)    = 0;
                    a_i(Model.ClusteredData.T < 5)      = 0;
                    gamma(Model.ClusteredData.T < 5)    = 0;
                    delta(Model.ClusteredData.T < 5)    = 0;
                end
            end        
    end

    % Unit Investment Annuity 
    C_i = r / (1 - (1 / (1 + r)^LT)) * InputData.EconomicData.DAC.ModuleInvestmentCost;
    % Unit Maintenence Cost
    C_m = InputData.EconomicData.DAC.MaintenanceCost * C_i;

    %% ------ DECISION VARIABLES
    Model.Ntot        = intvar(1);                      % Total number of modules
    Model.output      = sdpvar(D, K);                   % CO2 Output
    Model.E_el_DAC    = sdpvar(D, K);                   % Minimal electric energy requirement
    Model.E_th_DAC    = sdpvar(D, K);                   % Minimal thermal energy requirement
    Model.input{1}    = sdpvar(D, K);                   % Electric energy requirement
    Model.input{2}    = sdpvar(D, K);                   % Thermal energy requirement
    Model.E_elspl     = sdpvar(D, K);                   % Electric split
    Model.E_DACtot    = sdpvar(D, K);                   % Total energy requirement of DAC
    Model.Non         = sdpvar(D, K);                   % Number of units turned on

    if ConstantOP
        Model.x           = binvar(nrOS, 1);            % Operational parameter combination
        Model.Non_tild    = intvar(D, K, nrOS);         % Auxiliary Variable (Non)
        Model.Ntot_tild   = intvar(D, K, nrOS);         % Auxiliary Variable (Ntot)
    else
        Model.z_tild      = sdpvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (Energy_tot)
        Model.Non_tild    = intvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (Non)
        Model.y           = binvar(D, K, nrPieces);         % Piece identifier for input-output relation
        Model.v           = binvar(D, K, nrPieces);         % Piece identifier for input-input relation
        Model.y_tild      = binvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (y)
        % Identifier for OS
        if Settings.HourlySwitching == 1 || RHadapt
            Model.x           = intvar(D, K, nrOS);
        else
            Model.x           = intvar(nrOS,1);
        end

        Model.w_tild      = sdpvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (Energy_tot, el)
        Model.Non_hat     = intvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (Non)
        Model.v_tild      = binvar(D, K, nrOS, nrPieces);   % Piece identifier for input-input_el relation
        Model.Ntot_hat    = intvar(D, K, nrPieces);         % Auxiliary Variable (Ntot)

        Model.Ntot_tild   = intvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (Ntot)
        Model.Ntot_hat    = intvar(D, K, nrOS, nrPieces);   % Auxiliary Variable (Ntot)
    end
    
    %% ------ CONSTRAINTS
    if ConstantOP
        % CO2 Output
        Model.constraints = [Model.output == sum(Model.Non_tild .* alpha,3)];
        
        % Energy Input
        Model.constraints = [Model.constraints, Model.E_el_DAC == sum(Model.Non_tild .* epsilon_el,3)];
        Model.constraints = [Model.constraints, Model.E_th_DAC == sum(Model.Non_tild .* epsilon_th,3)];
        Model.constraints = [Model.constraints, Model.E_DACtot == Model.E_el_DAC +  Model.E_th_DAC];
        
        % Calculate Non
        Model.constraints = [Model.constraints, Model.Non == sum(Model.Non_tild,3)];
        
        % Sum constraint on x
        Model.constraints = [Model.constraints, sum(Model.x) == 1];
        
        % N_tot needs to be larger zero
        Model.constraints = [Model.constraints, Model.Ntot >= 0];
        
        for i = 1:nrOS
            %First constraint
            Model.constraints = [Model.constraints, 0 <= Model.Non_tild(:,:,i) <= Model.Ntot_tild(:,:,i)];
            %Second constraint
            Model.constraints = [Model.constraints, Model.size_min * Model.x(i) <= Model.Ntot_tild(:,:,i) <= Model.size_max * Model.x(i)];
            %Third constraint
            Model.constraints = [Model.constraints, Model.Ntot - Model.size_max*(1-Model.x(i)) <= Model.Ntot_tild(:,:,i) <= Model.Ntot];
        end
    else
        % CO2 Output
        Model.constraints = [Model.output == sum(sum(alpha .* Model.z_tild + beta .* Model.Non_tild,4),3)];
        % Electrical energy demand
        Model.constraints = [Model.constraints, Model.E_el_DAC == sum(sum(gamma .* Model.w_tild + delta .* Model.Non_hat,4),3)];

        % Calculate total energy input
        Model.constraints = [Model.constraints, Model.E_DACtot == sum(sum(Model.z_tild,4),3)];
        Model.constraints = [Model.constraints, Model.E_DACtot == sum(sum(Model.w_tild,4),3)];

        % Allow only one active piece in PW function
        Model.constraints = [Model.constraints, sum(sum(Model.y_tild,4),3) == 1];
        Model.constraints = [Model.constraints, sum(sum(Model.v_tild,4),3) == 1];

        % PWL for input-output relation
        Model.constraints = [Model.constraints, b_i_1 .* Model.Non_tild <= Model.z_tild];
        Model.constraints = [Model.constraints, Model.z_tild <= b_i .* Model.Non_tild];

        % PWL for input-input relation
        Model.constraints = [Model.constraints, a_i_1 .* Model.Non_hat <= Model.w_tild];
        Model.constraints = [Model.constraints, Model.w_tild <= a_i .* Model.Non_hat];

        % Size additional constraints
        Model.constraints = [Model.constraints, 0 <= Model.Non_tild];
        Model.constraints = [Model.constraints, Model.Non_tild <= Model.Ntot_tild];
        Model.constraints = [Model.constraints, 0 <= Model.Non_hat];
        Model.constraints = [Model.constraints, Model.Non_hat <= Model.Ntot_hat];

        Model.constraints = [Model.constraints, Model.size_min * Model.y_tild <= Model.Ntot_tild];
        Model.constraints = [Model.constraints, Model.Ntot_tild <= Model.size_max * Model.y_tild];
        Model.constraints = [Model.constraints, Model.size_min * Model.v_tild <= Model.Ntot_hat];
        Model.constraints = [Model.constraints, Model.Ntot_hat <= Model.size_max * Model.v_tild];


        for os = 1:nrOS
            if Settings.HourlySwitching == 1 || RHadapt
                for pie = 1:nrPieces
                    Model.constraints = [Model.constraints, 0 <= Model.y_tild(:, :, os, pie) <= Model.x(:,:,os)];
                    Model.constraints = [Model.constraints, 0 <= Model.v_tild(:, :, os, pie) <= Model.x(:,:,os)];
                    Model.constraints = [Model.constraints, Model.y(:, :, pie) - (1 - Model.x(:,:,os)) <= Model.y_tild(:, :, os, pie) <= Model.y(:, :, pie)];
                    Model.constraints = [Model.constraints, Model.v(:, :, pie) - (1 - Model.x(:,:,os)) <= Model.v_tild(:, :, os, pie) <= Model.v(:, :, pie)];

                    Model.constraints = [Model.constraints, Model.Ntot - Model.size_max * (1 - Model.y_tild(:, :, os, pie)) <= Model.Ntot_tild(:, :, os, pie) <= Model.Ntot];
                    Model.constraints = [Model.constraints, Model.Ntot - Model.size_max * (1 - Model.v_tild(:, :, os, pie)) <= Model.Ntot_hat(:, :, os, pie) <= Model.Ntot];
                    
                    Model.constraints = [Model.constraints, Model.Non - Model.Ntot + Model.Ntot_tild(:, :, os, pie) <= Model.Non_tild(:, :, os, pie)];
                    Model.constraints = [Model.constraints, Model.Non_tild(:, :, os, pie) <= Model.Non];
                    Model.constraints = [Model.constraints, Model.Non - Model.Ntot + Model.Ntot_hat(:, :, os, pie) <= Model.Non_hat(:, :, os, pie)];
                    Model.constraints = [Model.constraints, Model.Non_hat(:, :, os, pie) <= Model.Non];
                end
            else
                for pie = 1:nrPieces
                    Model.constraints = [Model.constraints, 0 <= Model.y_tild(:, :, os, pie) <= Model.x(os)];
                    Model.constraints = [Model.constraints, 0 <= Model.v_tild(:, :, os, pie) <= Model.x(os)];
                    Model.constraints = [Model.constraints, Model.y(:, :, pie) - (1 - Model.x(os)) <= Model.y_tild(:, :, os, pie) <= Model.y(:, :, pie)];
                    Model.constraints = [Model.constraints, Model.v(:, :, pie) - (1 - Model.x(os)) <= Model.v_tild(:, :, os, pie) <= Model.v(:, :, pie)];

                    Model.constraints = [Model.constraints, Model.Ntot - Model.size_max * (1 - Model.y_tild(:, :, os, pie)) <= Model.Ntot_tild(:, :, os, pie) <= Model.Ntot];
                    Model.constraints = [Model.constraints, Model.Ntot - Model.size_max * (1 - Model.v_tild(:, :, os, pie)) <= Model.Ntot_hat(:, :, os, pie) <= Model.Ntot];
                    
                    Model.constraints = [Model.constraints, Model.Non - Model.Ntot + Model.Ntot_tild(:, :, os, pie) <= Model.Non_tild(:, :, os, pie)];
                    Model.constraints = [Model.constraints, Model.Non_tild(:, :, os, pie) <= Model.Non];
                    Model.constraints = [Model.constraints, Model.Non - Model.Ntot + Model.Ntot_hat(:, :, os, pie) <= Model.Non_hat(:, :, os, pie)];
                    Model.constraints = [Model.constraints, Model.Non_hat(:, :, os, pie) <= Model.Non];
                end
            end

        end

        % Electric and thermal energy demand
        Model.constraints = [Model.constraints, Model.E_th_DAC == Model.E_DACtot - Model.E_el_DAC];
    end
    
    if TradeOff==1
        Model.constraints = [Model.constraints, Model.input{1} == Model.E_el_DAC + Model.E_elspl];
        Model.constraints = [Model.constraints, Model.input{2} == Model.E_DACtot - Model.E_el_DAC - eta_elth * Model.E_elspl];
        Model.constraints = [Model.constraints, Model.input{1} >= 0];
        Model.constraints = [Model.constraints, Model.input{2} >= 0];
        Model.constraints = [Model.constraints, Model.E_elspl >= 0];
    else
        Model.constraints = [Model.constraints, Model.input{1} == Model.E_el_DAC];
        Model.constraints = [Model.constraints, Model.input{2} == Model.E_th_DAC];
        Model.constraints = [Model.constraints, Model.E_elspl == 0];
    end
    
    % CO2 Production constraint
    if Settings.ConstantDemand == 1
        Model.constraints = [Model.constraints,...
            Model.output >= Model.ClusteredData.Demand];
    else
        Model.constraints = [Model.constraints,...
            sum(sum(Model.output(Model.ClusteredData.idx,:))) == sum(sum(Model.ClusteredData.Demand(Model.ClusteredData.idx,:)))];
    end
    %----- DEFINE OBJECTIVE
    switch Settings.OptimizationType
        case 1
        Model.objective = ...
            sum(sum(Model.input{1}(Model.ClusteredData.idx,:) .* Model.ClusteredData.P_el(Model.ClusteredData.idx,:))) ...
            + sum(sum(Model.input{2}(Model.ClusteredData.idx,:) .* Model.ClusteredData.P_th(Model.ClusteredData.idx,:))) ...
            + Model.Ntot * C_i ... 
            + Model.Ntot * C_m ;
        case 2
            Model.objective = ...
            sum(sum(Model.input{1}(Model.ClusteredData.idx,:) .* Model.ClusteredData.CO2FactorEl(Model.ClusteredData.idx,:))) ...
            + sum(sum(Model.input{2}(Model.ClusteredData.idx,:) .* Model.ClusteredData.CO2FactorTh(Model.ClusteredData.idx,:)));
        case 3
            Model.objective = ...
            sum(sum(Model.input{1}(Model.ClusteredData.idx,:) .* Model.ClusteredData.P_el(Model.ClusteredData.idx,:))) ...
            + sum(sum(Model.input{2}(Model.ClusteredData.idx,:) .* Model.ClusteredData.P_th(Model.ClusteredData.idx,:))) ...
            + Model.Ntot * C_i ... 
            + Model.Ntot * C_m ;

            Model.constraints = [Model.constraints,...
                Settings.EmissionConstraint >= ...
                sum(sum(Model.input{1}(Model.ClusteredData.idx,:) .* Model.ClusteredData.CO2FactorEl(Model.ClusteredData.idx,:)))...
              + sum(sum(Model.input{2}(Model.ClusteredData.idx,:) .* Model.ClusteredData.CO2FactorTh(Model.ClusteredData.idx,:)))...
            ];
    end
    if nargin == 8
        Model.constraints = [Model.constraints,...
                    Model.objective <= upperlimit];
    end
            
    solveroutput = optimize(Model.constraints,Model.objective,Solver_Options);
    

    %----- WRITE OUTPUT
    ModelResults.Solveroutput.solvertime = solveroutput.solvertime;
    ModelResults.Solveroutput.info = solveroutput.info;
    ModelResults.Solveroutput.problem = solveroutput.problem;
%     ModelResults.Solveroutput.mipgap = solveroutput.solveroutput.result.mipgap;
    ModelResults.Model.RH  = RH;
    ModelResults.Model.T   = T;
    ModelResults.Model.P_el  = P_el;
    ModelResults.Model.P_th  = P_th;
    ModelResults.Model.ClusteredData = Model.ClusteredData;
    ModelResults.TotalCosts = sum(sum(value(Model.input{1}(Model.ClusteredData.idx,:)) .* Model.ClusteredData.P_el(Model.ClusteredData.idx,:))) ...
            + sum(sum(value(Model.input{2}(Model.ClusteredData.idx,:)) .* Model.ClusteredData.P_th(Model.ClusteredData.idx,:))) ...
            + value(Model.Ntot) * C_i ... 
            + value(Model.Ntot) * C_m;
    ModelResults.TotalEmissions = sum(sum(value(Model.input{1}(Model.ClusteredData.idx,:)) .* Model.ClusteredData.CO2FactorEl(Model.ClusteredData.idx,:))) ...
            + sum(sum(value(Model.input{2}(Model.ClusteredData.idx,:)) .* Model.ClusteredData.CO2FactorTh(Model.ClusteredData.idx,:)));
    ModelResults.DACSize = value(Model.Ntot);
    ModelResults.DemandEl = sum(sum(value(Model.input{1}(Model.ClusteredData.idx,:))));
    ModelResults.DemandTh = sum(sum(value(Model.input{2}(Model.ClusteredData.idx,:))));
    ModelResults.EmissionsCaptured.Total = sum(sum(value(Model.output(Model.ClusteredData.idx,:))));
    ModelResults.Costs.Total = ModelResults.TotalCosts;
    ModelResults.Costs.Investment = value(Model.Ntot * C_i);
    ModelResults.Costs.Maintanance = value(Model.Ntot * C_m);
    ModelResults.Costs.Electricity = sum(sum(value(Model.input{1}(Model.ClusteredData.idx,:)) .* Model.ClusteredData.P_el(Model.ClusteredData.idx,:)));
    ModelResults.Costs.Heat = sum(sum(value(Model.input{2}(Model.ClusteredData.idx,:)) .* Model.ClusteredData.P_th(Model.ClusteredData.idx,:)));
    
    if SaveDetails
        ModelResults.Model.CO2FactorEl = CO2FactorEl;
        ModelResults.Model.CO2FactorTh = CO2FactorTh;
        ModelResults.Model.SolverOptions = Solver_Options;
        ModelResults.Model.Settings = Settings;
        ModelResults.EmissionsCaptured.TypDays = value(Model.output);
        ModelResults.EmissionsCaptured.FullRes = value(Model.output(Model.ClusteredData.idx,:));
        ModelResults.Input.Electricity.TypDays = value(Model.input{1});
        ModelResults.Input.Electricity.FullRes = value(Model.input{1}(Model.ClusteredData.idx,:));
        ModelResults.Input.Heat.TypDays = value(Model.input{2});
        ModelResults.Input.Heat.FullRes = value(Model.input{2}(Model.ClusteredData.idx,:));
        ModelResults.Output.TypDays = value(Model.output);
        ModelResults.Output.FullRes = value(Model.output(Model.ClusteredData.idx,:));
        if ConstantOP
            ModelResults.Operation.N_on = value(sum(Model.Non_tild,3));
            ModelResults.Operation.x = value(Model.x)' * linspace(1,size(value(Model.x),1),size(value(Model.x),1))';
            OP = Input_DAC(Input_DAC.OSIdentifier1 == ModelResults.Operation.x,:);
            OP = OP(1,:);
            temp.adTime = OP.adTime;
            temp.prodTime = OP.prodTime;
            temp.purgeTime = OP.purgeTime;
            temp.cycleTime = OP.cycleTime;
            temp.pressure = OP.pressure;
            temp.prodTemp = OP.prodTemp;
            temp.purgeDeltaTemp = OP.purgeDeltaTemp;
            temp.airVolume = OP.airVolume;
        else
            ModelResults.Operation.E_DACtot.TypDays = value(Model.E_DACtot);
            ModelResults.Operation.E_DACtot.FullRes = value(Model.E_DACtot(Model.ClusteredData.idx,:));
            ModelResults.Operation.y.TypDays = value(Model.y);
            if Settings.HourlySwitching == 1 || RHadapt >0
                ModelResults.Operation.x.TypDays = value(Model.x);
            else
                ModelResults.Operation.x.TypDays = value(Model.x)' *linspace(1,size(value(Model.x),1),size(value(Model.x),1))';
            end
                ModelResults.Operation.N_on = value(sum(sum(Model.Non_tild,3),4));
                x = Input_DAC.RH;
                y = Input_DAC.T;
                z = Input_DAC.E_tot;

                % adsorption time
                v = Input_DAC.adTime;
                F = scatteredInterpolant(x,y,z,v);
                T_d  = ModelResults.Model.ClusteredData.T;
                T_d(T_d<=min(Input_DAC.T)) = min(Input_DAC.T);
                T_d(T_d>=max(Input_DAC.T)) = max(Input_DAC.T);
                temp.adTime = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % production time
                v = Input_DAC.prodTime;
                F = scatteredInterpolant(x,y,z,v);
                temp.prodTime = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % purge time
                v = Input_DAC.purgeTime;
                F = scatteredInterpolant(x,y,z,v);
                temp.purgeTime = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % cycle time
                v = Input_DAC.cycleTime;
                F = scatteredInterpolant(x,y,z,v);
                temp.cycleTime = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % pressure
                v = Input_DAC.pressure;
                F = scatteredInterpolant(x,y,z,v);
                temp.pressure = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % production temperature
                v = Input_DAC.prodTemp;
                F = scatteredInterpolant(x,y,z,v);
                temp.prodTemp = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % purge temperature
                v = Input_DAC.purgeDeltaTemp;
                F = scatteredInterpolant(x,y,z,v);
                temp.purgeDeltaTemp = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);

                % air volume
                v = Input_DAC.airVolume;
                F = scatteredInterpolant(x,y,z,v);
                temp.airVolume = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on);
                ModelResults.Operation.OperationalParameters = temp;
                
                % water production
                v = Input_DAC.Water_Out;
                F = scatteredInterpolant(x,y,z,v);
                temp.waterProduction = F(ModelResults.Model.ClusteredData.RH, T_d, ModelResults.Operation.E_DACtot.TypDays ./ ModelResults.Operation.N_on) .* ModelResults.Operation.N_on;
                ModelResults.Operation.OperationalParameters = temp;
                
                if RHadapt >0
                    ModelResults.Operation.RH_choice = sum(value(Model.x) .* RH_choices, 3);
                    T_choices(isnan(T_choices)) = 0;
                    ModelResults.Operation.T_choice = sum(value(Model.x) .* T_choices, 3);
                    Water_cons(isnan(Water_cons)) = 0;
                    Water_cons = sum(value(Model.x) .* Water_cons, 3);
                    ModelResults.Operation.WaterConsumption = temp.airVolume .* Water_cons * 3600 * 3.4 * 10^6;
                    ModelResults.Operation.NetWaterProduction = -ModelResults.Operation.WaterConsumption + temp.waterProduction;
                    ModelResults.Operation.IncreasedRH = (ModelResults.Operation.RH_choice ~= ModelResults.Model.ClusteredData.RH);
                end
        end
        
    end
    toc
end
