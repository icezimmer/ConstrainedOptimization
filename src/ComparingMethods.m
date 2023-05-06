function [table_results, table_solutions] = ComparingMethods(Q, q, P, date, off_the_shelves, eps_R, max_steps)
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

if nargin < 5 % no eps, max_steps
    eps_R = 1e-5;
    max_steps = 1000;
elseif nargin == 5 % no max_steps
    max_steps = 1000;
end

f_star = Optimum(Q, q, P);

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
Histories = cell(1, 3);

% Away-step Frank-Wolfe type algorithm
[x_min, f_min, elapsed_time, num_steps, type, step_size_method, converging, feasible, duality_gap, history] = AwayStepFrankWolfe(Q, q, P, eps_R, max_steps, false, false, date, f_star);
Method = cat(1, Method, type);
Step_Size = cat(1, Step_Size, step_size_method);
Minimum = cat(1, Minimum, f_min);
Time = cat(1, Time, elapsed_time);
Steps = cat(1, Steps, num_steps);
Converging = cat(1, Converging, converging);
Feasible = cat(1, Feasible, feasible);
Duality_Gap = cat(1, Duality_Gap, duality_gap);
Solutions = cat(2, Solutions, x_min);
Histories{1} = history;

% Quadratic Programming algorithms
for i = 1:length(off_the_shelves)
    [x_min, f_min, elapsed_time, num_steps, type, step_size_method, converging, feasible, duality_gap, history] = OffTheShelf(Q, q, P, off_the_shelves(i), max_steps, eps_R, f_star);
    Method = cat(1, Method, type);
    Step_Size = cat(1, Step_Size, step_size_method);
    Minimum = cat(1, Minimum, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
    Histories{1+i} = history;
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

label_column = {'FW_Exact', 'Direct'};
label_column = cat(2,label_column,string(off_the_shelves));
table_solutions = array2table(Solutions, 'VariableNames', label_column);

colors = ['k','b','r','g','y'];
iterative_algorithms = cat(2,Method{1},off_the_shelves);
gcf = figure('Name', 'Comparison');
hold on
for i = 1:length(iterative_algorithms)
    gap = abs(Histories{i} - f_star) / max(1,abs(f_star));
    semilogy(0:length(gap)-1, gap, 'Color', colors(i),'DisplayName', iterative_algorithms(i));
end
hold off
set(gca, 'YScale', 'log');
title('Comparison')
xlabel('step')
ylabel('error')
legend('Location','best')
saveas(gcf, fullfile('results', date, 'comparison.png'))

end