function [table_results, table_solutions] = ComparingMethods(Q, q, P, date, eps_R, max_steps)
%{
Grid search on line search method and momentum coefficient.
Input:
    Q                  : (matrix) nxn positive semi-definite
    q                  : (vector) of length n
    P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    date               : (string) date for saving the results
    eps_R              : (float) max relative error
    max_steps          : (integer) stop criterion for Frank Wolfe (max number of steps for Frank Wolfe)
Output:
    table_results   : (table) results obtained by the methods. They are sorted by the minimum, duality_gap, time
        and number of steps
    table_solutions : (table) table of solutions (x_min) for all the methods
%}

if nargin < 5 % no eps, max_steps
    eps_R = 0.1;
    max_steps = 1000;
elseif nargin == 5 % no max_steps
    max_steps = 1000;
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
step_size_methods = "Exact";
[x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, step_size_methods, eps_R, max_steps, false, false, date);
Method = cat(1, Method, type);
Step_Size = cat(1, Step_Size, step_size_method);
Minimum = cat(1, Minimum, f_min);
Time = cat(1, Time, elapsed_time);
Steps = cat(1, Steps, num_steps);
Converging = cat(1, Converging, converging);
Feasible = cat(1, Feasible, feasible);
Duality_Gap = cat(1, Duality_Gap, duality_gap);
Solutions = cat(2, Solutions, x_min);

% Quadratic Programming algorithms
algorithms = ["interior-point-convex", "active-set"];
for i = 1:length(algorithms)
    [x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = QuadraticProgramming(Q, q, P, algorithms(i), max_steps);
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

table_solutions = array2table(Solutions, 'VariableNames', {'FW_Exact', 'Direct', 'interior-point-convex', 'active-set'});

end