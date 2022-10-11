%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and number of simplices
n = 50; K = ceil(0.2*n);
% Kernel dimension and spectral radius of the matrix Q (it must be > 0)
dim_ker = 0; spectral_radius = 1;
% Density of the matrix Q
density = 1;
% Minimum value, maximum value and number of zero in the vector q
min_q = -5; max_q = 5; zero_q = 0;
% Seed for the random generator
seed = 4;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, K, dim_ker, spectral_radius, density, min_q, max_q, zero_q, seed);

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, date)

% Save the variables (not necessary)
SaveVariables(Q, q, P, date)

% Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
eps_R = 1e-8; max_steps = 1e2;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, date, eps_R, max_steps);

% Save the results
SaveResults(table_results, table_solutions, date)