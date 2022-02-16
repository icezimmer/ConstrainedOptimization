function [table_results, table_solutions] = ComparingMethods(Q, q, P, x_start, eps, max_steps, eps_ls, tomography, optimization_curve, date)
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
    tomography         : (logical) plot or not the tomography for each step
    optimization_curve : (logical) plot or not the optimization curve
Output:
    table_results   : (table) results obtained by the methods. They are sorted by the minimum, duality_gap, time
        and number of steps
    table_solutions : (table) table of solutions (x_min) for all the methods
%}

[n, ~] = size(Q);

step_size_methods = ["Default", "LBM", "QBM", "NM"];

Method = zeros(0,1);
Minimum = zeros(0,1);
Time = zeros(0,1);
Step_Size = zeros(0,1);
Steps = zeros(0,1);
Converging = zeros(0,1);
Feasible = zeros(0,1);
Duality_Gap = zeros(0,1);

Solutions = zeros(n, 0);

for i = 1:length(step_size_methods)
    [x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, step_size_methods(i), tomography, optimization_curve, date);
    Method = [Method; type];
    Step_Size = [Step_Size; step_size_method];
    Minimum = [Minimum; f_min];
    Time = [Time; elapsed_time];
    Steps = [Steps; num_steps];
    Converging = [Converging; converging];
    Feasible = [Feasible; feasible];
    Duality_Gap = [Duality_Gap; duality_gap];
    Solutions = [Solutions, x_min];
end

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


[x_min, f_min, elapsed_time, type, step_size_method, num_steps, converging, feasible, duality_gap] = QuadraticProgramming(Q, q, P, x_start, max_steps);
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

table_solutions = array2table(Solutions, 'VariableNames', {'FW_Default', 'FW_LBM', 'FW_QBM', 'FW_NM', 'Direct', 'QuadProg'});

end