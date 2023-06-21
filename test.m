%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension, number of simplices, force or not
% non-point-simplices, and the fraction of active constraints respect the
% solution
n = 1000; K = 500; force_non_point_simplices = true; actv = 0;
% Kernel dimension, spectral radius of the matrix Q, and minimum
% eigenvalue (considered only if dim_ker>0)
dim_ker = 0; spectral_radius = 1; lambda_min = 1;
% Density of the matrix Q
density = 1;
% Norm of the vector q
% norm_q = 1000;

% Seed for the random generator
seed = 7;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, num_vertex, date] = GenerateInstance(n, K, force_non_point_simplices, actv, dim_ker, spectral_radius, lambda_min, density, seed);
%norm(q)

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, num_vertex, date)

% Save the variables
SaveVariables(Q, q, P, date)

% Stopping criteria for the Frank Wolfe method: max relative error and max number of steps for Frank Wolfe
eps_R = 1e-10; max_steps = 1e4;
% Define the step size selection method: "Away-step" or "Standard"
variant = "Away-step";

% Plot or not the tomography for each iteration
tomography = false;
% Plot or not the error curve
error_plot = true;
% Perform the Frank-Wolfe algorithm
[x_min, f_min, elapsed_time, num_steps, method, variant, converging, feasible, duality_gap, history] = FrankWolfe(Q, q, P, variant, eps_R, max_steps, tomography, error_plot, date);
