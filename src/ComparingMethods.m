function [table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, varargin)
%{
Comparason between the FW algorithm and two off-the-shelf methods (interior-point-convex and active-set)
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

numvarargs = length(varargin);
if numvarargs > 5
    error('myfuns:ComparingMethods:TooManyInputs', ...
        'requires at most 5 optional inputs');
end

% set defaults for optional inputs
optargs = {false,1e-6,1e-8,1e-10,1000};
optargs(1:numvarargs) = varargin;
[dual_comparison, eps_RDG, eps_RT, eps_RE, max_steps] = optargs{:};

% Pre-compute the optimum by the oracle
[~, f_star] = Oracle(Q, q, P);

[n, ~] = size(Q);

Method = zeros(0,1);
Minima = zeros(0,1);
Time = zeros(0,1);
Variants = zeros(0,1);
Steps = zeros(0,1);
Error = zeros(0,1);
Converging = zeros(0,1);
Feasible = zeros(0,1);
Duality_Gap = zeros(0,1);

Solutions = zeros(n, 0);
Histories = cell(1, 0);

% Frank-Wolfe type algorithms
for i = 1:length(frank_wolfe_variants)
    [x_min, f_min, elapsed_time, num_steps, type, variant, err, converging, feasible, duality_gap, history] = FrankWolfe(Q, q, P, frank_wolfe_variants(i), eps_RDG, eps_RE, max_steps, false, false, date, f_star);
    Method = cat(1, Method, type);
    Variants = cat(1, Variants, variant);
    Minima = cat(1, Minima, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Error = cat(1, Error, err);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
    Histories = cat(2, Histories, history);
end

% Quadratic Programming algorithms
for i = 1:length(off_the_shelves)
    [x_min, f_min, elapsed_time, num_steps, type, variant, err, converging, feasible, duality_gap, history] = QuadraticProgramming(Q, q, P, off_the_shelves(i), max_steps, eps_RT, eps_RE, f_star);
    Method = cat(1, Method, type);
    Variants = cat(1, Variants, variant);
    Minima = cat(1, Minima, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Error = cat(1, Error, err);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
    Histories = cat(2, Histories, history);
end

table_results = table(Method, Variants, Minima, Time, Steps, Error, Converging, Feasible, Duality_Gap);
table_results = sortrows(table_results, {'Error', 'Time', 'Steps'});

label_column = string(frank_wolfe_variants);
label_column = cat(2,label_column,string(off_the_shelves));
table_solutions = array2table(Solutions, 'VariableNames', label_column);

colors = ["#0072BD","#D95319","#77AC30","#EDB120","#7E2F8E"];
iterative_algorithms = cat(2,append('FW-',string(frank_wolfe_variants)),off_the_shelves);
PlotBenchmark(Histories, f_star, iterative_algorithms, dual_comparison, colors, date)

end