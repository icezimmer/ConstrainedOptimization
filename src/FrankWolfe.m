function [x_min, f_min, elapsed_time, num_steps, converging, duality_gap] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, step_size_method, tomography, optimization_curve)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain. 
Input:
    Q                  : (matrix) nxn positive semi-definite
    q                  : (vector) of length n
    P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start            : (vector) start point
    eps                : (float) stop criterion (max duality_gap for Frank Wolfe)
    max_steps          : (integer) stop criterion (max number of steps for Frank Wolfe)
    eps_ls             : (float) stop criterion for the line search
    step_size_method   : (string) method for line search
    tomography         : (logical) plot or not the tomography for each step
    optimization_curve : (logical) plot or not the optimization curve
Output:
    x_min        : (vector) argmin of the function
    f_min        : (vector) min of the function
    elapsed_time : (float) time elapsed for the computation
    num_steps    : (integer) number of steps for the convergence
    converging   : (string) method converge or not
    duality_gap  : (float) opposite value of the scalar product between
        the descent direction and the gradient
%}

if isequal(step_size_method,'NM')
    disp('LineSearch by Newton Method')
elseif isequal(step_size_method,'LBM')
    disp('LineSearch by Linear Bisection Method')
elseif isequal(step_size_method,'QBM')
    disp('LineSearch by Quadratic Bisection Method')
elseif isequal(step_size_method,'Default')
    disp('Default Step Size Selection')
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

if optimization_curve
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
    if optimization_curve
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

% Check the convergence of the algorithm
if(duality_gap >  eps || ~Domain(x_min, P))
    converging = "No";
else
    converging = "Yes";
end

% Final duality_gap (duality gap) 
duality_gap = E(end);

if tomography
    disp(['it. ', num2str(i), ', f(x) = ', num2str(f_min)])
end

% Plot the optimization curve
if optimization_curve
    fx(i+1) = f_min;
    PlotOptimizationCurve(fx, E, step_size_method)
end

end