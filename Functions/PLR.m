function [x_bp, y_bp, SSR, alpha, beta] = PLR(X, Y, nr_bp, opt_fixed_xbp)
%PLR Piecewise linear regressions with floating breakpoints.
%   [x_bp, y_bp, SSR, alpha, beta] = PLR(X, Y, nr_bp) finds the n = nr_bp
%   breakpoints that minimize the sum of squared residuals between the data
%   and the fitted linear pieces. The fitted function is assumed to be
%   continuous. x_bp and y_bp return the x and y coordinates of the
%   breakpoints. SSR is the resultsing sum of squared residuals, and alpha
%   and beta respective slope and intercept parameters of the pieces. If
%   the input argument opt_fixed_xbp is populated, the function assumes
%   fixed x-values for the breakpoints
    
    if (nargin <= 3) 
        % Set initial bp vectors
        x_bp0 = linspace(min(X), max(X), nr_bp);
        y_bp0 = linspace(min(Y), max(Y), nr_bp);
        free_bp0 = [x_bp0(2:end-1) y_bp0];
    
        % Minimize the SSR
        [free_bp, SSR] = fminsearch(@(free_bp) SSR_func(free_bp, X, Y, nr_bp), free_bp0);
    
        % Calculate breakpoints
        x_bp = [min(X), free_bp(1:nr_bp-2), max(X)];
        y_bp = free_bp(nr_bp-1:end);
    else
        % Set initial bp vectors
        y_bp0 = linspace(min(Y), max(Y), nr_bp);
        free_bp0 = y_bp0;
    
        % Minimize the SSR
        options = optimset('TolFun',1e-3);
        [free_bp, SSR] = fminsearch(@(free_bp) SSR_func(free_bp, X, Y, nr_bp, opt_fixed_xbp), free_bp0, options);
    
        % Calculate breakpoints
        x_bp = opt_fixed_xbp;
        y_bp = free_bp;
    end

    % Calculate alpha and beta
    alpha = zeros(nr_bp-1,1);
    beta  = zeros(nr_bp-1,1);
    for i=1:nr_bp-1
        alpha(i) = (y_bp(i+1)-y_bp(i))/(x_bp(i+1)-x_bp(i));
        beta(i)  = y_bp(i)-alpha(i)*x_bp(i);
    end

end

function SSR = SSR_func(free_bp, x_data, y_data, nr_bp, opt_fixed_xbp)
%SSR_func Calculating the SSR of a piecewise defined linear function.
%   SSR = SSR_func(free_bp, x_data, y_data, nr_bp) returns the sum of
%   squared residuals with the x values of the first and last breakpoint
%   being fixed a the max and min of the data. free_bp is a vector of
%   remaining breakpoint coordinates (i.e. all y values and remaining x
%   values. nr_bp specifies the number of breakpoints including the corner
%   points
    
    % All breakpoint variables
    if (nargin <= 4) 
        x_bp = [min(x_data), free_bp(1:nr_bp-2), max(x_data)];
        y_bp = free_bp(nr_bp-1:end);
    else
        x_bp = opt_fixed_xbp;
        y_bp = free_bp;
    end
    
    
    % Calculate alpha and beta
    alpha = zeros(nr_bp-1,1);
    beta = zeros(nr_bp-1,1);
    for i=1:nr_bp-1
        alpha(i) = (y_bp(i+1)-y_bp(i))/(x_bp(i+1)-x_bp(i));
        beta(i)  = y_bp(i)-alpha(i)*x_bp(i);
    end
    
    % Calculate identifier for each datapoint per piece
    for i=1:nr_bp-1
        if i==nr_bp-1
            identifier(:,i) = (x_bp(i) <= x_data & x_data <= x_bp(i+1));
        else
            identifier(:,i) = (x_bp(i) <= x_data & x_data < x_bp(i+1));
        end
    end
    
    % Calculate Residuals
    for i=1:nr_bp-1
        y_hat(:,i) = alpha(i)*x_data+beta(i);
    end
    
    % Calculate SSR
    SSR = sum((sum(y_hat .* identifier,2)-y_data) .^2);
end
