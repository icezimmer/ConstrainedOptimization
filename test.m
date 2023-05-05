%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and number of simplices
n = 1000; K = ceil(0.2*n);
% Kernel dimension and spectral radius of the matrix Q (it must be > 0)
dim_ker = 100; spectral_radius = 10;
% Density of the matrix Q
density = 1;
% Minimum value, maximum value and number of zero in the vector q
min_q = -5; max_q = 5; zero_q = 0;
% Seed for the random generator
seed = 7;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, K, dim_ker, spectral_radius, density, min_q, max_q, zero_q, seed);

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, date)

% Save the variables
SaveVariables(Q, q, P, date)

% Stopping criteria for the Frank Wolfe method: max relative error and max number of steps for Frank Wolfe
eps_R = 1e-6; max_steps = 2*1e3;
% Define the step size selection method: 'Exact' or 'Standard'
step_size_method = 'Exact';

% Plot or not the tomography for each iteration
tomography = false;
% Plot or not the error curve
error_plot = true;
% Perform the Frank-Wolfe algorithm
[x_min, f_min, elapsed_time, num_steps, method, step_size_method, converging, feasible, duality_gap] = Prova(Q, q, P, step_size_method, eps_R, max_steps, tomography, error_plot, date);
