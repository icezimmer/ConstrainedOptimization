function [table_results, table_solutions] = ComparingMethods(Q, q, P, date, x_start, eps, max_steps, eps_ls)
%{
Grid search on line search method and momentum coefficient.
Input:
    Q                  : (matrix) nxn positive semi-definite
    q                  : (vector) of length n
    P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    date               : (string) date for saving the results
    x_start            : (vector) start point
    eps                : (float) stop criterion for Frank Wolfe (max Duality_gap for Frank Wolfe)
    max_steps          : (integer) stop criterion for Frank Wolfe (max number of steps for Frank Wolfe)
    eps_ls             : (float) stop criterion for the line search 
Output:
    table_results   : (table) results obtained by the methods. They are sorted by the minimum, duality_gap, time
        and number of steps
    table_solutions : (table) table of solutions (x_min) for all the methods
%}

if nargin < 5 % no x_start, eps, max_steps, eps_ls
    x_start = zeros(length(Q), 1);
    [K, ~] = size(P);
    for k = 1 : K
        i = find(P(k,:), 1);
        x_start(i) = 1;
    end
    eps = 0.1;
    max_steps = 1000;
    eps_ls = 0.01;
elseif nargin == 5 % no eps, max_steps, eps_ls
    eps = 0.1;
    max_steps = 1000;
    eps_ls = 0.01;
elseif nargin == 6 % no max_steps, eps_ls
    max_steps = 1000;
    eps_ls = 0.01;
elseif nargin == 7 % no eps_ls
    eps_ls = 0.01;
end

[n, ~] = size(Q);

Method = zeros(0,1);
Minimum = zeros(0,1);
Time = zeros(0,1);
Step_Size = zeros(0,1);
Steps = zeros(0,1);
Converging = zeros(0,1);
Feasible = zeros(0,1);
Duality_Gap = zeros(0,1);

Solutions = zeros(n, 0);

% Frank-Wolfe type algorithms
step_size_methods = ["Default", "LBM", "QBM", "NM"];
for i = 1:length(step_size_methods)
    [x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, date, x_start, eps, max_steps, eps_ls, step_size_methods(i));
    Method = cat(1, Method, type);
    Step_Size = cat(1, Step_Size, step_size_method);
    Minimum = cat(1, Minimum, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
end

% Quadratic Programming algorithms
algorithms = ["interior-point-convex", "active-set"];
for i = 1:length(algorithms)
    [x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = QuadraticProgramming(Q, q, P, algorithms(i), x_start, max_steps);
    Method = cat(1, Method, type);
    Step_Size = cat(1, Step_Size, step_size_method);
    Minimum = cat(1, Minimum, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
end

% Direct Computation
[x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = DirectComputation(Q, q, P);
Method = [Method; type];
Step_Size = [Step_Size; step_size_method];
Minimum = [Minimum; f_min];
Time = [Time; elapsed_time];
Steps = [Steps; num_steps];
Converging = [Converging; converging];
Feasible = [Feasible; feasible];
Duality_Gap = [Duality_Gap; duality_gap];
Solutions = [Solutions, x_min];


table_results = table(Method, Step_Size, Minimum, Time, Steps, Converging, Feasible, Duality_Gap);
table_results = sortrows(table_results, {'Minimum', 'Duality_Gap', 'Time', 'Steps'});

table_solutions = array2table(Solutions, 'VariableNames', {'FW_Default', 'FW_LBM', 'FW_QBM', 'FW_NM', 'Direct', 'interior-point-convex', 'active-set'});

end