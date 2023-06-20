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
if numvarargs > 2
    error('myfuns:ComparingMethods:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {1e-5,1000};
optargs(1:numvarargs) = varargin;
[eps_R, max_steps] = optargs{:};

f_star = Optimum(Q, q, P);

[n, ~] = size(Q);

Method = zeros(0,1);
Minima = zeros(0,1);
Time = zeros(0,1);
Variants = zeros(0,1);
Steps = zeros(0,1);
Converging = zeros(0,1);
Feasible = zeros(0,1);
Duality_Gap = zeros(0,1);

Solutions = zeros(n, 0);
Histories = cell(1, 0);

% Away-step Frank-Wolfe type algorithm
for i = 1:length(frank_wolfe_variants)
    [x_min, f_min, elapsed_time, num_steps, type, variant, converging, feasible, duality_gap, history] = FrankWolfe(Q, q, P, frank_wolfe_variants(i), eps_R, max_steps, false, false, date, f_star);
    Method = cat(1, Method, type);
    Variants = cat(1, Variants, variant);
    Minima = cat(1, Minima, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
    Histories = cat(2, Histories, history);
end

% Quadratic Programming algorithms
for i = 1:length(off_the_shelves)
    [x_min, f_min, elapsed_time, num_steps, type, variant, converging, feasible, duality_gap, history] = OffTheShelf(Q, q, P, off_the_shelves(i), max_steps, eps_R, f_star);
    Method = cat(1, Method, type);
    Variants = cat(1, Variants, variant);
    Minima = cat(1, Minima, f_min);
    Time = cat(1, Time, elapsed_time);
    Steps = cat(1, Steps, num_steps);
    Converging = cat(1, Converging, converging);
    Feasible = cat(1, Feasible, feasible);
    Duality_Gap = cat(1, Duality_Gap, duality_gap);
    Solutions = cat(2, Solutions, x_min);
    Histories = cat(2, Histories, history);
end

table_results = table(Method, Variants, Minima, Time, Steps, Converging, Feasible, Duality_Gap);
table_results = sortrows(table_results, {'Minima', 'Duality_Gap', 'Time', 'Steps'});

label_column = string(frank_wolfe_variants);
label_column = cat(2,label_column,string(off_the_shelves));
table_solutions = array2table(Solutions, 'VariableNames', label_column);

colors = ['k','b','r','g','y'];
iterative_algorithms = cat(2,append('FW-',string(frank_wolfe_variants)),off_the_shelves);
gcf = figure('Name', 'Comparison');
hold on
for i = 1:length(iterative_algorithms)
    gap = abs(Histories{i} - f_star) / max(1,abs(f_star));
%     if isequal(iterative_algorithms(i), 'FW-Away-step')
%         r = 1 / sqrt(2); % max rate of linear convergence for AFW
%         semilogy(0:length(gap)-1, gap(1) * r.^(0:length(gap)-1), 'Color', 'm', 'DisplayName', "Max-Rate-AFW");
%     end
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