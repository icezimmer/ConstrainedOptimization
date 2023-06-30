%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain (Minkowslìki sum of unitary simplices).
%}

addpath src

% Space dimension, number of simplices, force or not non-point-simplices, and the fraction of active constraints respect the solution
n = 100; K = 10; force_non_point_simplices = true; actv = 0.5;
% Kernel dimension, spectral radius of the matrix Q (considered only if dim_ker<n), and minimum eigenvalue (considered only if dim_ker=0)
dim_ker = 10; spectral_radius = 2; lambda_min = 1;
% Density of the matrix Q
density = 1;

% Seed for the random generator
seed = 1;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, num_vertex, date] = GenerateInstance(n, K, force_non_point_simplices, actv, dim_ker, spectral_radius, lambda_min, density, seed);

% Save the parameters
SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_ker, spectral_radius, lambda_min, density, date)

% Save the matrices
SaveMatrices(Q, q, P, date)

% Stopping criteria for the Frank Wolfe method: max relative error and max number of steps for Frank Wolfe
eps_R = 1e-10; max_steps = 1e5;
% Define the step size selection method: "Away-step" or "Standard"
variant = "Away-step";

% Plot or not the tomography for each iteration
tomography = false;
% Plot or not the error curve
error_plot = true;
% Perform the Frank-Wolfe algorithm
[x_min, f_min, elapsed_time, num_steps, method, variant, converging, feasible, duality_gap, history] = FrankWolfe(Q, q, P, variant, eps_R, max_steps, tomography, error_plot, date);

% Save the results
SaveTestResults(x_min, f_min, elapsed_time, num_steps, method, variant, converging, feasible, duality_gap, history, date)