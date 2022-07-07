function [x_min, f_min, elapsed_time, method,  step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, date, x_start, eps, max_steps, eps_ls, step_size_method, tomography, optimization_curve, convergence_rate)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain. 
Input:
    Q                  : (matrix) nxn positive semi-definite
    q                  : (vector) of length n
    P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    date               : (string) date for saving the results
    x_start            : (vector) starting point
    eps                : (float) stop criterion (max duality_gap for Frank Wolfe)
    max_steps          : (integer) stop criterion (max number of steps for Frank Wolfe)
    eps_ls             : (float) stop criterion for the line search
    step_size_method   : (string) method for line search
    tomography         : (logical) plot or not the tomography for each step
    optimization_curve : (logical) plot or not the optimization curve
    convergence_rate   : (logical) plot or not the log-log optimization curve
Output:
    x_min        : (vector) argmin of the function
    f_min        : (vector) min of the function
    elapsed_time : (float) time elapsed for the computation
    num_steps    : (integer) number of steps for the convergence
    converging   : (logical) method converge or not
    feasible     : (logical) solution is feasible or not
    duality_gap  : (float) opposite value of the scalar product between
        the descent direction and the gradient
%}

if nargin < 5 % no x_start, eps, max_steps, eps_ls, step_size_method, tomography, optimization_curve, convergence_rate
    x_start = zeros(length(Q), 1);
    [K, ~] = size(P);
    for k = 1 : K
        i = find(P(k,:), 1);
        x_start(i) = 1;
    end
    eps = 0.1;
    max_steps = 1000;
    eps_ls = 0.01;
    step_size_method = 'Default';
    tomography = false;
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 5 % no eps, max_steps, eps_ls, step_size_method, tomography, optimization_curve, convergence_rate
    eps = 0.1;
    max_steps = 1000;
    eps_ls = 0.01;
    step_size_method = 'Default';
    tomography = false;
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 6 % no max_steps, eps_ls, step_size_method, tomography, optimization_curve, convergence_rate
    max_steps = 1000;
    eps_ls = 0.01;
    step_size_method = 'Default';
    tomography = false;
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 7 % no eps_ls, step_size_method, tomography, optimization_curve, convergence_rate
    eps_ls = 0.01;
    step_size_method = 'Default';
    tomography = false;
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 8 % no step_size_method, tomography, optimization_curve, convergence_rate
    step_size_method = 'Default';
    tomography = false;
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 9 % no tomography, optimization_curve, convergence_rate
    tomography = false;
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 10 % no optimization_curve, convergence_rate
    optimization_curve = true;
    convergence_rate = false;
elseif nargin == 11 % no convergence_rate
    convergence_rate = false;
end

method = "FW";

if isequal(step_size_method,'NM')
    disp('Frank-Wolfe algorithm with LineSearch by Newton Method')
elseif isequal(step_size_method,'LBM')
    disp('Frank-Wolfe algorithm with LineSearch by Linear Bisection Method')
elseif isequal(step_size_method,'QBM')
    disp('Frank-Wolfe algorithm with LineSearch by Quadratic Bisection Method')
elseif isequal(step_size_method,'Default')
    disp('Frank-Wolfe algorithm with Default Step Size Selection')
end

% Force q and x to be column vectors
q = q(:);
x = x_start(:);

% Function f
f = @(x) x'*Q*x + q'*x;

tic
% First iteration
i = 0;

[d, ~, duality_gap] = LinearApproximationMinimizer(Q, q, P, x);

% Line search
[alpha, alpha_start] = StepSizeSelection(Q, q, x, d, eps_ls, i, step_size_method);

if tomography
    Tomography(Q, q, x, d, step_size_method, alpha, alpha_start, i, duality_gap)
end

if (optimization_curve || convergence_rate)
    fx(i+1) = f(x);
    E(i+1) = duality_gap;
end

% Upgrade the vector x
x = NewPoint(x, d, alpha);

% New iteration
i = i + 1;

% Iterate until convergence
while (duality_gap >  eps && i < max_steps)
    [d, ~, duality_gap] = LinearApproximationMinimizer(Q, q, P, x);
    
    [alpha, alpha_start] = StepSizeSelection(Q, q, x, d, eps_ls, i, step_size_method);
    
    % Plot tomography
    if tomography
        Tomography(Q, q, x, d, step_size_method, alpha, alpha_start, i, duality_gap)
    end
    
    % Append new value of the funtion and the duality_gap for the optimization curve
    if (optimization_curve || convergence_rate)
        fx(i+1) = f(x);
        E(i+1) = duality_gap;
    end
    
    % Upgrade the point
    x = NewPoint(x, d, alpha);
    
    i = i + 1;
end

% Elapsed time for computation
elapsed_time = toc;

% Minimum belonging to the domain
x_min = x;
% Minimum value of the function
f_min = f(x);

% Number of steps
num_steps = i;

% Convergence of the algorithm
converging = (duality_gap <=  eps);

% Feasibility of the solution
feasible = CheckDomain(x_min, P);

if tomography
    disp(['it. ', num2str(i), ', f(x) = ', num2str(f_min)])
end

% Append the value of the function (last iteration)
if (optimization_curve || convergence_rate)
    fx(i+1) = f_min;
end

% Plot the optimization curve
if optimization_curve
    PlotOptimizationCurve(fx, E, step_size_method, date)
end

% Plot the log-log optimization curve
if convergence_rate
    PlotConvergenceRate(fx, E, step_size_method, date)
end    

end