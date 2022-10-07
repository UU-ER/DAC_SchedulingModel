%% Load Data
% Load sorbent performance data
load('DAC_Data_5_CT80.mat')
% Load Optimization Information
load('Settings.mat')
% Load Input data
load('Input_DAC.mat');
% Load Solver options
load('Solver_Options.mat')
% You can add external solvers like Gurobi or CPLEX as a solver here
solver_options.solver = '';


%% Adapt Settings
% Cluster data to D typical days
Settings.D = 100;
% Turn on/off (1/0) ohmic heating
Settings.ElThTradeOff = 1;
% Assume flat CO2 production profile or not (0/1)
Settings.ConstantDemand = 0;
% Assume constant operating parameters or not (0/1)
Settings.ConstantOP = 0;
% Water spraying of inlet air on/off (1/0)
Settings.RHadapt = 0;
% If water spraying is on, how many discrete points?
Settings.RHdiscret = 11;
% Minimize cost (1), emissions (2), or cost at emission constraint (3)
Settings.OptimizationType = 1;
% Set emission constraint
Settings.EmissionConstraint = 0;

%% Adapt Input Data
% Change the Input struct respectively

%% Set RH and T vectors
RH = ones(1,8760)*50;
T = ones(1,8760)*15;

%% Run optimization
ones(1,8760)*50;
results = DAC_Standalone_6(solver_options, Settings, Input, RH, T, DAC_Data, 1);
