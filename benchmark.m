%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and kernel dimension of the matrix Q
n = 10000; dim_ker = 3;
% Spectral radius of the matrix Q (it must be > 0)
spectral_radius = 10;
% Density of the matrix Q
dns = 0.0001;
% Minimum value, maximum value and number of zero in the vector q
min_q = -5; max_q = 5; zero_q = 0;
% Seed for the random generator
seed = 0;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, x_start, K_plus, date] = GenerateInstance(n, seed, dim_ker, spectral_radius, dns, min_q, max_q, zero_q);
%[Q, q, P, x_start, K_plus, date] = GenerateInstance2(n, seed, 10, 3, 0.1, min_q, max_q, zero_q);

% Save the parameters
PrintParameters(n, dim_ker, spectral_radius, K_plus, date)

% Save the variables (not necessary)
%PrintVariables(Q, q, P, date)

% Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
eps = 0.1; max_steps = 100;
% Stop criterion for the line search methods
eps_ls = 0.01;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, x_start, eps, max_steps, eps_ls, date);

% Save the results
PrintResults(table_results, table_solutions, date)

