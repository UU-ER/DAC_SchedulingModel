function [r_b_i_1, r_b_i, r_alpha, r_beta, r_a_i_1, r_a_i, r_gamma, r_delta, r_b_iy_1, r_b_iy, r_a_iy_1, r_a_iy] = DACOutput_4(T, RH, nr_bp, DAC_Data)
%DACOutput_4 generates the performance parameters of a DAC module based on
%the data provided in DAC_Data as a linearily approximated function
%   T and RH are temperature and relative humidity at which the performance parameters are
%   evaluated. The function then
%   (1) linearly interpolates the DAC_Data for the given relative humidity
%   (RH) and temperature (T) resulting in two (or more) pareto lines for
%   two different operational strategies (at the time of coding, it is
%   minimizing energy and minimizing exergy).
%   (2) To further simplify the subsequent MILP problem, the interpolated 
%   data is reduced to a linearily approximated function with n = nr_bp 
%   breakpoints (see PLR function, n includes the endpoints).
%   Resulting parameters are:
%
%   r_b_i_1 and r_b_i are the interval borders of the respective piece in
%   the piecewiese defined function for the Total Energy Input-CO2 output 
%   relation.
%   r_alpha and r_beta are respectivly slope and intercept parameters of the
%   respective piecewise defined function.
%   r_a_i_1 and r_a_i are the interval borders of the respective piece in
%   the piecewiese defined function for the Total Energy Input-electricity  
%   input relation.
%   r_gamma and r_delta are respectivly slope and intercept parameters of the
%   respective piecewise defined function.
%
% TODO: 
% - Implement water production
        
if size(RH,1) ~= size(T,1) || size(RH,2) ~= size(T,2)
   error('Size of temperature and humidity matrix are not equal.')
end
disp('DAC module: Starting to derive pareto lines');
N = size(RH,1);   % number of days per year 
K = size(RH,2);   % number of daily intervals

% Number of pareto points
nrPar = max(DAC_Data.Point);
% Number of operating strategies
nrOS  = unique(DAC_Data.OS);
nrOS  = size(nrOS,1);
if nrOS > 1
   error('You can only calculate the performance parameters for one OS')
end

% Create interpolated pareto points
Out    = zeros(N,K,1,nrPar);
E_tot  = zeros(N,K,1,nrPar);
E_el   = zeros(N,K,1,nrPar);
OS     = zeros(N,K,1,nrPar);
Point  = zeros(N,K,1,nrPar);

% Create Resulting matrices
r_b_i_1    = zeros(N,K,nr_bp-1);
r_b_i      = zeros(N,K,nr_bp-1);
r_alpha    = zeros(N,K,nr_bp-1);
r_beta     = zeros(N,K,nr_bp-1);
r_a_i_1    = zeros(N,K,nr_bp-1);
r_a_i      = zeros(N,K,nr_bp-1);
r_gamma    = zeros(N,K,nr_bp-1);
r_delta    = zeros(N,K,nr_bp-1);

% Prevent the extrapolation of temperature data
T(T<=min(DAC_Data.T)) = min(DAC_Data.T);
T(T>=max(DAC_Data.T)) = max(DAC_Data.T);

% Make all interpolated paretos
% Loop through all pareto points
for pt = 1:nrPar
    data = DAC_Data(DAC_Data.Point == pt,:);
    x = data.RH;
    y = data.T;
    % CO2 Output
    v = data.CO2_Out;
    F = scatteredInterpolant(x,y,v);
    Out(:,:,1,pt)       = F(RH, T);
    % Total energy requirements
    v = data.E_tot;
    F = scatteredInterpolant(x,y,v);
    E_tot(:,:,1,pt)     = F(RH, T);
    % Electrical energy requirements
    v = data.E_el;
    F = scatteredInterpolant(x,y,v);
    E_el(:,:,1,pt)      = F(RH, T);
    % Point identifier
    Point(:,:,1,pt)     = pt;
end

% Loop through each timeslot and lineary approximate the interpolated
% pareto lines of the last step
for n=1:N
    for k=1:K
        if ~isnan(T(n,k))
            % PLR for Input-Output relation
            X = squeeze(E_tot(n,k,1,:));
            Y = squeeze(Out(n,k,1,:));
            [x_bp_os, y_bp_os, ~, ~, ~] = PLR(X, Y, nr_bp);
            y_bp_os(y_bp_os<=1e-6)=0;
            alpha_os = zeros(nr_bp-1,1);
            beta_os  = zeros(nr_bp-1,1);
            for i=1:nr_bp-1
                if (x_bp_os(i+1)-x_bp_os(i)) == 0
                    alpha_os(i,1) = 0;
                    beta_os(i,1) = 0;
                else
                    alpha_os(i,1) = (y_bp_os(i+1)-y_bp_os(i))/(x_bp_os(i+1)-x_bp_os(i));
                    beta_os(i,1)  = y_bp_os(i)-alpha_os(i)*x_bp_os(i);
                end
            end
            alpha(:,1) = alpha_os;
            beta(:,1)  = beta_os;
            for piece = 1:nr_bp-1
                b_i_1(piece,1) = x_bp_os(piece);
                b_i(piece,1) = x_bp_os(piece+1);
                b_iy_1(piece,1) = y_bp_os(piece);
                b_iy(piece,1) = y_bp_os(piece+1);
            end
            % PLR for Electricity-Total Energy relation
            X = squeeze(E_tot(n,k,1,:));
            Y = squeeze(E_el(n,k,1,:));
            [x_bp_os, y_bp_os, ~, ~, ~] = PLR(X, Y, nr_bp);
            x_bp_os(x_bp_os<=1e-6)=0;
            y_bp_os(y_bp_os<=1e-6)=0;
            x_bp_os(y_bp_os<=1e-6)=0;
            y_bp_os(x_bp_os<=1e-6)=0;
            gamma_os = zeros(nr_bp-1,1);
            delta_os  = zeros(nr_bp-1,1);
            for i=1:nr_bp-1
                if (x_bp_os(i+1)-x_bp_os(i)) == 0
                    gamma_os(i,1) = 0;
                    delta_os(i,1) = 0;
                else
                    gamma_os(i,1) = (y_bp_os(i+1)-y_bp_os(i))/(x_bp_os(i+1)-x_bp_os(i));
                    delta_os(i,1)  = y_bp_os(i)-gamma_os(i)*x_bp_os(i);
                end
            end
            gamma(:,1) = gamma_os;
            delta(:,1)  = delta_os;
            for piece = 1:nr_bp-1
                a_i_1(piece,1) = x_bp_os(piece);
                a_i(piece,1) = x_bp_os(piece+1);
                a_iy_1(piece,1) = y_bp_os(piece);
                a_iy(piece,1) = y_bp_os(piece+1);
            end
            % Write results to resulting matrices 
            r_b_i_1(n,k,:)    = reshape(b_i_1,[],1);
            r_b_i(n,k,:)      = reshape(b_i,[],1);
            r_b_iy_1(n,k,:)   = reshape(b_iy_1,[],1);
            r_b_iy(n,k,:)     = reshape(b_iy,[],1);
            r_alpha(n,k,:)    = reshape(alpha,[],1);
            r_beta(n,k,:)     = reshape(beta,[],1);
            r_a_i_1(n,k,:)    = reshape(a_i_1,[],1);
            r_a_i(n,k,:)      = reshape(a_i,[],1);
            r_a_iy_1(n,k,:)   = reshape(a_iy_1,[],1);
            r_a_iy(n,k,:)     = reshape(a_iy,[],1);
            r_gamma(n,k,:)    = reshape(gamma,[],1);
            r_delta(n,k,:)    = reshape(delta,[],1);
        else
            r_b_i_1(n,k,:)    = 0;
            r_b_i(n,k,:)      = 0;
            r_b_iy_1(n,k,:)   = 0;
            r_b_iy(n,k,:)     = 0;
            r_alpha(n,k,:)    = 0;
            r_beta(n,k,:)     = 0;
            r_a_i_1(n,k,:)    = 0;
            r_a_i(n,k,:)      = 0;
            r_a_iy_1(n,k,:)   = 0;
            r_a_iy(n,k,:)     = 0;
            r_gamma(n,k,:)    = 0;
            r_delta(n,k,:)    = 0;
        end
    end
end
disp('DAC module: Input Data processed')
end
