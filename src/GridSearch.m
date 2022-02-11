function [x_min, best, grid_search] = GridSearch(Q, q, P, x_start, eps, max_steps, eps_ls, candidates, tomography, optimization_curve)
%{
Grid search on line search method and momentum coefficient.
Input:
    Q                  : (matrix) nxn positive semi-definite
    q                  : (vector) of length n
    P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start            : (vector) start point
    eps                : (float) stop criterion for Frank Wolfe (max Duality_gap for Frank Wolfe)
    max_steps          : (integer) stop criterion for Frank Wolfe (max number of steps for Frank Wolfe)
    eps_ls             : (float) stop criterion for the line search
    candidates         : (Map object) container for wich keys are the hyperparameters, and values the possible values of each parameter 
    tomography         : (logical) plot or not the tomography for each step
    optimization_curve : (logical) plot or not the optimization curve
Output:
    x_min       : (vector) solution in the domain
    best        : (table) best method and parameters
    grid_search : (table) all methods sorted
%}

step_size_methods = candidates('ls_method');
step_size_methods = step_size_methods(:);

Minimum = zeros(0,1);
Time = zeros(0,1);
Steps = zeros(0,1);
Step_Size = step_size_methods;
Converging = zeros(0,1);
Duality_Gap = zeros(0,1);
Solution = zeros(0, length(x_start));

for i = 1:length(step_size_methods)
    [x_min, f_min, elapsed_time, num_steps, converging, duality_gap] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, step_size_methods(i), tomography, optimization_curve);
    Minimum = [Minimum; f_min];
    Time = [Time; elapsed_time];
    Steps = [Steps; num_steps];
    Converging = [Converging; converging];
    Duality_Gap = [Duality_Gap; duality_gap];
    Solution = [Solution; x_min'];
end

grid_search = table(Minimum, Time, Steps, Step_Size, Converging, Duality_Gap, Solution);
grid_search = sortrows(grid_search, {'Minimum', 'Duality_Gap', 'Time', 'Steps'});

x_min = grid_search(1, end);
x_min = table2array(x_min);
x_min = x_min';

grid_search = grid_search(:, 1 : end-1);

best = grid_search(1, :);

end