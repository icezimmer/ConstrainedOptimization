%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Seed for the random generator
seed = 0;
% Load the sparse matrix Q and generate the sparse vector q and the starting point x_start
[Q, q, P, x_start, date] = GenerateFromExternalData('bcsstk09.mtx', seed);


% Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
eps = 0.1; max_steps = 99;
% Stop criterion for the line search method
eps_ls = 0.01;
% Define the step size selection method
step_size_method = 'Default';

% Plot or not the tomography for each iteration
tomography = false;
% Plot or not the optimization curve
optimization_curve = true;
% Plot or not the log-log otimization curve (for the convergence rate)
convergence_rate = true;
[x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, step_size_method, tomography, optimization_curve, convergence_rate, date);