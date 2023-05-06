function [x_min, f_min, elapsed_time, num_steps, method,  step_size_method, converging, feasible, duality_gap, history] = AwayStepFrankWolfe(Q, q, P, eps_R, max_steps, tomography, error_plot, date, f_star)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain. 
Input:
    Q                : (matrix) nxn positive semi-definite
    q                : (vector) of length n
    P                : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    step_size_method : (string) method for line search
    eps_R            : (float) maximum relative error for stop condition
    max_steps        : (integer) stop criterion (max number of steps for Frank Wolfe)
    tomography       : (logical) plot or not the tomography for each step
    error_plot       : (logical) plot or not the error curve
    date             : (string) date for saving the results
    f_star           : (float) optimal value of the function in the domain
Output:
    x_min            : (vector) argmin of the function
    f_min            : (vector) min of the function
    elapsed_time     : (float) time elapsed for the computation
    method           : (string)
    step_size_method : (string)
    num_steps        : (integer) number of steps for the convergence
    converging       : (logical) method converge or not
    feasible         : (logical) solution is feasible or not
    duality_gap      : (float) opposite value of the scalar product between
        the descent direction and the gradient
%}

if nargin < 4
    eps_R = 1e-5;
    max_steps = 10000;
    tomography = false;
    error_plot = false;
    date = string(datetime('now','TimeZone','local','Format','d-MMM-y_HH:mm:ss'));
    mkdir(fullfile('results',date));
    f_star = Optimum(Q, q(:), P);
elseif nargin == 4
    max_steps = 10000;
    tomography = false;
    error_plot = false;
    date = string(datetime('now','TimeZone','local','Format','d-MMM-y_HH:mm:ss'));
    mkdir(fullfile('results',date));
    f_star = Optimum(Q, q(:), P);
elseif nargin == 5
    tomography = false;
    error_plot = false;
    date = string(datetime('now','TimeZone','local','Format','d-MMM-y_HH:mm:ss'));
    mkdir(fullfile('results',date));
    f_star = Optimum(Q, q(:), P);
elseif nargin == 6
    error_plot = false;
    date = string(datetime('now','TimeZone','local','Format','d-MMM-y_HH:mm:ss'));
    mkdir(fullfile('results',date));
    f_star = Optimum(Q, q(:), P);
elseif nargin == 7
    date = string(datetime('now','TimeZone','local','Format','d-MMM-y_HH:mm:ss'));
    mkdir(fullfile('results',date));
    f_star = Optimum(Q, q(:), P);
elseif nargin == 8
    f_star = Optimum(Q, q(:), P);
end

step_size_method = "Exact";
method = "AFW";

% Construct the starting point for the Frank-Wolfe algorithm
x_start = StartingPoint(P);

% Force q and x to be column vectors
q = q(:);
x = x_start(:);

% Function f
f = @(x) x'*Q*x + q'*x;

disp('Away-step Frank-Wolfe algorithm with Exact Line Search')

[indices, partition] = PartitionDomain(P);

tic
i = 0;
history = f(x);
duality_gap = NaN;
E = zeros(0,1);
% Iterate until convergence
while (~StoppingCriteria(history(end), f_star, eps_R) && i < max_steps)

    [d, ~, duality_gap, d_a, ~, duality_gap_a, alpha_max] = LinearizationMinimizerAFW(Q, q, x, indices, partition);
    
    if duality_gap >= duality_gap_a
        alpha_max = 1;
    else
        d = d_a;
        duality_gap = duality_gap_a;
    end

    %alpha_max

    alpha = StepSizeSelection(Q, d, duality_gap, alpha_max, i, step_size_method);
    
    % Plot tomography
    if tomography
        Tomography(Q, q, x, d, alpha, alpha_max, i, duality_gap)
    end
    
    % Append new value of the duality_gap for the error plot
    if error_plot
        E(i+1) = duality_gap;
    end

    % Upgrade the point
    x = NewPoint(x, d, alpha);

    % Append new value of the funtion
    history = cat(2, history, f(x));
    
    i = i + 1;
end

% Elapsed time for computation
elapsed_time = toc;

% Minimum belonging to the domain
x_min = x;
% Minimum value of the function
f_min = history(end);

% Number of steps
num_steps = i;

% Feasibility of the solution
feasible = CheckDomain(x_min, P);

% Convergence of the algorithm
converging = StoppingCriteria(f(x), f_star, eps_R) && feasible;

if tomography
    disp(['it. ', num2str(i), ', f(x) = ', num2str(f_min)])
end

% Plot the error
if error_plot
    PlotErrorCurve(history, f_star, E, step_size_method, date)
end  

end