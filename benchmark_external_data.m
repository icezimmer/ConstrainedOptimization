%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Seed for the random generator
seed = 0;
% Load the sparse matrix Q and generate the sparse vector q and the starting point x_start
[Q, q, P, x_start, minima, n, dim_ker, spectral_radius, K_plus, date] = GenerateFromExternalData('bcsstm13.mtx', seed);

% Save the parameters
PrintParameters(n, dim_ker, spectral_radius, K_plus, date)

% Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
eps = 0.1; max_steps = 500;
% Stop criterion for the line search methods
eps_ls = 0.01;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, x_start, eps, max_steps, eps_ls, date);

% Save the results
PrintResults(table_results, table_solutions, date)

