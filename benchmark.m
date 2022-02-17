%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and kernel dimension of the matrix Q
n = 500; dim_ker = 0;
% Minimum and maximum of the strictly positive eigenvalues of the matrix Q 
min_eig = 0.25; max_eig = 10;
% Minimum value, maximum value and number of zero in the vector q
min_q = -5; max_q = 5; zero_q = 0;
% Seed for the random generator
seed = 0;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, x_start, minima, spectral_radius, K_plus, date] = GenerateInstance(n, dim_ker, min_eig, max_eig, min_q, max_q, zero_q, seed);

% Save the parameters
PrintParameters(n, dim_ker, spectral_radius, K_plus, date)

% Save the variables
PrintVariables(Q, q, P, date)

% Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
eps = 0.1; max_steps = 1000;
% Stop criterion for the line search methods
eps_ls = 0.01;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, x_start, eps, max_steps, eps_ls, date);

% Save the results
PrintResults(table_results, table_solutions, date)

