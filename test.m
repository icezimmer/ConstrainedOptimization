%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and number of simplices
n = 1; K = 1; simplifyable = true;
% Kernel dimension and spectral radius of the matrix Q (it must be > 0)
dim_ker = 0; spectral_radius = 1000;
% Density of the matrix Q
density = 1;
% Minimum value, maximum value and number of zero in the vector q
min_q = 0; max_q = 10; zero_q = 0;
% Seed for the random generator
seed = 2;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, K, simplifyable, dim_ker, spectral_radius, density, min_q, max_q, zero_q, seed);

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, date)

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
