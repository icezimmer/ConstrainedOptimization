%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a constrained convex domain.
%}

addpath src

% Space dimension and kernel dimension of the matrix Q
n = 100; dim_ker = 10;
% Minimum and maximum of the strictly positive eigenvalues of the matrix Q 
min_eig = 1; max_eig = 10;
% Minimum value, maximum value and number of zero in the vector q
min_q = 3; max_q = 9; zero_q = 3;
% Seed for the random generator
seed = 0;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, x_start, minima] = GenerateInstance(n, dim_ker, min_eig, max_eig, min_q, max_q, zero_q, seed);

% Check if the global minimum exists and if it satisfies the constraints
is_in = false;
if minima
    [x_min, f_min] = GlobalMinimum(Q, q);
    is_in = Domain(x_min, P);
    if is_in
        x_min
        f_min
    end
end

% If the global minimum is not in the domain or it doesn't exist compute a grid search the Franke-Wolfe method
if ~is_in
    % Stoping criteria for the Frank Wolfe method: max error and max number of steps for Frank Wolfe
    eps = 0.1; max_steps = 500;
    % Stop criterion for the line search methods
    eps_ls = 0.01;
    
    % Define a Map object for the grid search (the line search methods and the beta values (momentum coefficients))
    candidates = containers.Map({'ls_method', 'beta'}, ...
        {["Default", "LBM", "QBM", "NM"], [0]});
    
    % Plot or not the tomography for each step and/or the optimization curve for each method
    tomography = false; curve = true;

    % Search the best solution
    [x_min, solution, grid_search] = GridSearch(Q, q, P, x_start, eps, max_steps, eps_ls, candidates, tomography, curve);
    
    % Print the solution x, the info about the best model and the grid search table
    x_min
    solution
    grid_search
end