function [x_min, best, grid_search] = GridSearch(Q, q, P, x_start, eps, max_steps, eps_ls, candidates, tomography, curve)
%{
Grid search on line search method and momentum coefficient.
Input:
    Q          : (matrix) nxn positive semi-definite
    q          : (vector) of length n
    P          : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start    : (vector) start point
    eps        : (float) stop criterion for Frank Wolfe (max Duality_gap for Frank Wolfe)
    max_steps  : (integer) stop criterion for Frank Wolfe (max number of steps for Frank Wolfe)
    eps_ls     : (float) stop criterion for the line search
    candidates : (Map object) container for wich keys are the hyperparameters, and values the possible values of each parameter 
    tomography : (logical) plot or not the tomography for each step
    curve      : (logical) plot or not the optimization curve
Output:
    x_min       : (vector) solution in the domain
    best        : (table) best method and parameters
    grid_search : (table) all methods sorted
%}

line_search_methods = candidates('ls_method');
line_search_methods = line_search_methods(:);

beta_values = candidates('beta');
beta_values = beta_values(:);

Minimum = zeros(0,1);
Time = zeros(0,1);
Steps = zeros(0,1);
LineSearch = repelem(line_search_methods, length(beta_values));
Momentum = repmat(beta_values, length(line_search_methods), 1);
Converging = zeros(0,1);
Duality_Gap = zeros(0,1);
Solution = zeros(0, length(x_start));

for i = 1:length(line_search_methods)
    for j = 1:length(beta_values)
        [x_min, f_min, elapsed_time, num_steps, converging, duality_gap] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, line_search_methods(i), beta_values(j), tomography, curve);
        Minimum = [Minimum; f_min];
        Time = [Time; elapsed_time];
        Steps = [Steps; num_steps];
        Converging = [Converging; converging];
        Duality_Gap = [Duality_Gap; duality_gap];
        Solution = [Solution; x_min'];
    end
end

grid_search = table(Minimum, Time, Steps, LineSearch, Momentum, Converging, Duality_Gap, Solution);
grid_search = sortrows(grid_search, {'Minimum', 'Duality_Gap', 'Time', 'Steps'});

x_min = grid_search(1, end);
x_min = table2array(x_min);
x_min = x_min';

grid_search = grid_search(:, 1 : end-1);

best = grid_search(1, :);

end