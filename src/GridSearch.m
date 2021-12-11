function [x_min, best, grid_search] = GridSearch(Q, q, P, x_start, eps, max_steps, eps_ls, line_search_methods, left, right, delta, tomography, curve)
%{
Grid search on line search method and momentum coefficient.
Input:
    Q          : (matrix) nxn positive semi-definite
    q          : (vector) of length n
    P          : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start    : (vector) start point
    eps        : (float) stop criterion for Frank Wolfe (max error for Frank Wolfe)
    max_steps  : (integer) stop criterion (max number of steps for Frank Wolfe)
    eps_ls     : (float) stop criterion for the line search
    left       : (string) minimum value for momentum coefficient
    right      : (float) maximum value for momentum coefficient 
    delta      : (float) with size of the sub intervals
    tomography : (logical) plot or not the tomography for each step
    curve       : (logical) plot or not the optimization curve
Output:
    x_min       : (vector) solution in the domain
    best        : (table) best method and parameters
    grid_search : (table) all methods sorted
%}

line_search_methods = line_search_methods(:);

if left == right
    beta_values = [left];
else
    beta_values = (left:delta:right)';
end
    
%line_search_methods = ["Trivial"; "LBM"; "QBM"; "NM"];

Minimum = zeros(0,1);
Time = zeros(0,1);
Steps = zeros(0,1);
LineSearch = repelem(line_search_methods, length(beta_values));
Momentum = repmat(beta_values, length(line_search_methods), 1);
Converging = zeros(0,1);
Error = zeros(0,1);
Solution = zeros(0, length(x_start));

for i = 1:length(line_search_methods)
    for j = 1:length(beta_values)
        [x_min, f_min, elapsed_time, num_steps, converging, error] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, line_search_methods(i), beta_values(j), tomography, curve);
        Minimum = [Minimum; f_min];
        Time = [Time; elapsed_time];
        Steps = [Steps; num_steps];
        Converging = [Converging; converging];
        Error = [Error; error];
        Solution = [Solution; x_min'];
    end
end

grid_search = table(Minimum, Time, Steps, LineSearch, Momentum, Converging, Error, Solution);
grid_search = sortrows(grid_search, {'Minimum', 'Time', 'Steps', 'Error'});

x_min = grid_search(1, end);
x_min = table2array(x_min);
x_min = x_min';

grid_search = grid_search(:, 1 : end-1);

best = grid_search(1, :);

end