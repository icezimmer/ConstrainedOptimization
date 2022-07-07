%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and kernel dimension of the matrix Q
n = 1000; dim_ker = 0;
% Spectral radius of the matrix Q (it must be > 0)
spectral_radius = 1000;
% Density of the matrix Q
density = 0.005;
% Minimum value, maximum value and number of zero in the vector q
min_q = -5; max_q = 5; zero_q = 0;
% Seed for the random generator
seed = 0;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, x_start, K_plus, date] = GenerateInstance(n, seed, dim_ker, spectral_radius, density, min_q, max_q, zero_q);

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, date)

% Save the variables
SaveVariables(Q, q, P, date)

% Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
eps = 0.1; max_steps = 1000;
% Stop criterion for the line search method
eps_ls = 0.01;
% Define the step size selection method
step_size_method = 'Default';

% Plot or not the tomography for each iteration
tomography = false;
% Plot or not the optimization curve
optimization_curve = false;
% Plot or not the log-log otimization curve (for the convergence rate)
convergence_rate = false;
% Perform the Frank-Wolfe algorithm
[x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, date, x_start, eps, max_steps, eps_ls, step_size_method, tomography, optimization_curve, convergence_rate);
