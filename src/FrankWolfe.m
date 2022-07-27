function [x_min, f_min, elapsed_time, method,  step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, step_size_method, eps_R, max_steps, tomography, error_plot, date)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain. 
Input:
    Q                  : (matrix) nxn positive semi-definite
    q                  : (vector) of length n
    P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    date               : (string) date for saving the results
    eps                : (float) stop criterion (max duality_gap for Frank Wolfe)
    max_steps          : (integer) stop criterion (max number of steps for Frank Wolfe)
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
if nargin < 4
    step_size_method = 'Default';
    eps_R = 1e-5;
    max_steps = 10000;
    tomography = false;
    error_plot = false;
    formatOut = 'mm_dd_yy_HH_MM_SS';
    date = datestr(now,formatOut);
    mkdir(fullfile('results',date));
elseif nargin == 4
    eps_R = 1e-5;
    max_steps = 10000;
    tomography = false;
    error_plot = false;
    formatOut = 'mm_dd_yy_HH_MM_SS';
    date = datestr(now,formatOut);
    mkdir(fullfile('results',date));
elseif nargin == 5
    max_steps = 10000;
    tomography = false;
    error_plot = false;
    formatOut = 'mm_dd_yy_HH_MM_SS';
    date = datestr(now,formatOut);
    mkdir(fullfile('results',date));
elseif nargin == 6
    tomography = false;
    error_plot = false;
    formatOut = 'mm_dd_yy_HH_MM_SS';
    date = datestr(now,formatOut);
    mkdir(fullfile('results',date));
elseif nargin == 7
    error_plot = false;
    formatOut = 'mm_dd_yy_HH_MM_SS';
    date = datestr(now,formatOut);
    mkdir(fullfile('results',date));
elseif nargin == 8
    formatOut = 'mm_dd_yy_HH_MM_SS';
    date = datestr(now,formatOut);
    mkdir(fullfile('results',date));
end

method = "FW";

if isequal(step_size_method,'Exact')
    disp('Frank-Wolfe algorithm with Exact Line Search')
elseif isequal(step_size_method,'Default')
    disp('Frank-Wolfe algorithm with Default Step Size Selection')
else
    error('Wrong Line Search name')
end

% Construct the starting point for the Frank-Wolfe algorithm
x_start = StartingPoint(P);

% Force q and x to be column vectors
q = q(:);
x = x_start(:);

% Function f
f = @(x) x'*Q*x + q'*x;

% Compute the minimum with off-the-shelf method
H = 2 * Q;
vec = q;
%Aeq = P_slim;
Aeq = P;
beq = ones(size(Aeq,1),1);
lb = zeros(size(Aeq,2),1);
options = optimoptions(@quadprog, 'Algorithm', "interior-point-convex", 'Display', 'off');
[~, f_star] = quadprog(H, vec, [], [], Aeq, beq, lb, [], [], options);

[indices, partition] = PartitionDomain(P);

tic
i = 0;
fx = zeros(0,1);
E = zeros(0,1);
% Iterate until convergence
%while (duality_gap >  eps && i < max_steps)
while (~StoppingCriteria(f(x), f_star, eps_R) && i < max_steps)

    [d, ~, duality_gap] = LinearApproximationMinimizer(Q, q, x, indices, partition);
    
    alpha = StepSizeSelection(Q, q, x, d, i, step_size_method);
    
    % Plot tomography
    if tomography
        Tomography(Q, q, x, d, alpha, i, duality_gap)
    end
    
    % Append new value of the funtion and the duality_gap for the error plot
    if error_plot
        fx(i+1) = f(x);
        E(i+1) = duality_gap;
    end

    % Upgrade the point
    x = NewPoint(x, d, alpha);
    
    i = i + 1;
end

% Elapsed time for computation
elapsed_time = toc;


% Enlarge the solution putting the fixed dimension
%x(indices) = x_slim;

% Minimum belonging to the domain
x_min = x;
% Minimum value of the function
f_min = f(x_min);

% Number of steps
num_steps = i;

% Convergence of the algorithm
converging = StoppingCriteria(f(x), f_star, eps_R);

% Feasibility of the solution
feasible = CheckDomain(x_min, P);

if tomography
    disp(['it. ', num2str(i), ', f(x) = ', num2str(f_min)])
end

% Append the value of the function (last iteration) and plot the error
if error_plot
    fx(i+1) = f(x);
    PlotOptimizationCurve(fx, f_star, E, step_size_method, date)
end  

end