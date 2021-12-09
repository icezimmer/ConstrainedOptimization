function [x_min, best, grid_search] = GridSearch(Q, q, P, x_start, eps, eps_ls, left, right, delta)
%{
Grid search on line search method and momentum coefficient.
Input:
    Q       : (matrix) nxn positive semi-definite
    q       : (vector) of length n
    P       : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start : (vector) start point
    eps     : (float) stop criterion for Frank Wolfe
    eps_ls  : (float) stop criterion for the line search
    left    : (string) minimum value for momentum coefficient
    right   : (float) maximum value for momentum coefficient 
    delta   : (float) with size of the sub intervals
Output:
    x_min       : (vector) solution in the domain
    best        : (table) best method and parameters
    grid_search : (table) all methods sorted
%}

if left == right
    beta_values = [left];
else
    beta_values = left:delta:right;
end
    
line_search_methods = ["LBM"; "QBM"; "NM"];

Minimum = zeros(0,1);
Time = zeros(0,1);
Steps = zeros(0,1);
LineSearch = repelem(line_search_methods, length(beta_values));
Momentum = repmat(beta_values, length(line_search_methods), 1);
Solution = zeros(0, length(x_start));

for i = 1:length(line_search_methods)
    for j = 1:length(beta_values)
        [x_min, f_min, elapsed_time, num_steps] = FrankWolfe(Q, q, P, x_start, eps, eps_ls, line_search_methods(i), beta_values(j));
        Minimum = [Minimum; f_min];
        Time = [Time; elapsed_time];
        Steps = [Steps; num_steps];
        Solution = [Solution; x_min'];
    end
end

grid_search = table(Minimum, Time, Steps, LineSearch, Momentum, Solution);
grid_search = sortrows(grid_search, {'Minimum', 'Time', 'Steps'});
x_min = grid_search(1, end);
x_min = table2array(x_min);
x_min = x_min';
grid_search = grid_search(:, 1 : end-1);
best = grid_search(1, :);