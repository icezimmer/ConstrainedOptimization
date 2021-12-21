%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a constrained convex domain.
%}

addpath src

% Space dimension and kernel dimension of the matrix Q
n = 80; dim_ker = 0;
% Minimum eigenvalue of the matrix Q and the maximum eigenvalue of the matrix Q
min_eig = 2; max_eig = 10;

% Minimum value, maximum value and number of zero in the vector q
min_q = 3; max_q = 9; zero_q = 3;

% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, x_start, minima] = Generate(n, dim_ker, min_eig, max_eig, min_q, max_q, zero_q);

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
    % Stop criteria: max error and max number of steps for Frank Wolfe
    eps = 0.1; max_steps = 1000;
    % Stop criterion for the line search methods
    eps_ls = 0.01;

    % List of line search methods to try ("Trivial", "LBM", "QBM", "NM")
    line_search_methods = ["QBM", "Trivial", "LBM", "NM"];
    % Minimum value of momentum coefficient, maximum value of momentum coefficient and width of the sub-intervals for the momentum coefficient
    left = 0; right = 0.1; delta = 0.01;

    % Plot or not the tomography for each step or the optimization curves for each step
    tomography = false; curve = false;

    % Search the best solution
    [x_min, solution, grid_search] = GridSearch(Q, q, P, x_start, eps, max_steps, eps_ls, line_search_methods, left, right, delta, tomography, curve);
    
    % Print the solution x, the info about the best model and the grid search table
    x_min
    solution
    grid_search
end