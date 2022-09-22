function [r_alpha, r_epsilon_el, r_epsilon_th, OS] = DACOutput_4_ConsPar(T, RH, DAC_Data)
%DACOutput_3 generates the performance parameters of a DAC module based on
%the data provided in DAC_Data as a linearily approximated function.
%Operating Parameters are assumed to be constant over the whole timeframe.
%   T and RH are temperature and relative humidity at which the performance parameters are
%   evaluated. The function then linearly interpolates the DAC_Data for the given relative humidity
%   (RH) and temperature (T).
%
%   r_alpha is the input-output parameter.
%   r_epsilon_el and r_epsilon_th are electric and thermal demand
%   parameters. 
        
if size(RH,1) ~= size(T,1) || size(RH,2) ~= size(T,2)
   error('Size of temperature and humidity matrix are not equal.')
end
disp('DAC module: Starting to derive pareto lines');
N = size(RH,1);   % number of days per year 
K = size(RH,2);   % number of daily intervals

% Number of operating strategies
nrOS  = unique(DAC_Data.OSIdentifier1);
nrOS  = size(nrOS,1);

% Create interpolated pareto points
Out    = zeros(N,K,nrOS);
E_el   = zeros(N,K,nrOS);
E_th   = zeros(N,K,nrOS);
OS     = zeros(N,K,nrOS);

% Create Resulting matrices
r_alpha      = zeros(N,K,nrOS);
r_epsilon_el = zeros(N,K,nrOS);
r_epsilon_th = zeros(N,K,nrOS);

% Prevent the extrapolation of temperature data
T(T<=min(DAC_Data.T)) = min(DAC_Data.T);
T(T>=max(DAC_Data.T)) = max(DAC_Data.T);

% Loop through all OP combinations

for OPC = 1:nrOS
    data = DAC_Data(DAC_Data.OSIdentifier1 == OPC,:);
    if size(unique(data.RH),1)==1 || size(unique(data.T),1)==1 || size(data,1) == 0
        Out(:,:,OPC) = zeros(N,K);
        E_el(:,:,OPC) = zeros(N,K);
        E_th(:,:,OPC) = zeros(N,K);
    else
        x = data.RH;
        y = data.T;
        % CO2 Output
        v = data.CO2_Out;
        F = scatteredInterpolant(x,y,v,'linear','none');
        Out(:,:,OPC)       = F(RH, T);
        % Electrical energy requirements
        v = data.E_el;
        F = scatteredInterpolant(x,y,v,'linear','none');
        E_el(:,:,OPC)      = F(RH, T);
        % Thermal energy requirements
        v = data.E_th;
        F = scatteredInterpolant(x,y,v,'linear','none');
        E_th(:,:,OPC)      = F(RH, T);
        % OS identifier
        OS(:,:,OPC)     = OPC;
    end
    
end
% Set all extrapolated values to zero (NaN -> 0)
Out(isnan(Out))=0;
E_el(isnan(E_el))=0;
E_th(isnan(E_th))=0;

% Derive alpha, epsilon
r_alpha = Out;
r_epsilon_el = E_el;
r_epsilon_th = E_th;

disp('DAC module: Input Data processed')
end
