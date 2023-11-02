function [table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, varargin)
    %{
        Comparason between FW algorithms and off-the-shelf methods (from quadprog)
        Input:
            Q                    : (matrix) nxn positive semi-definite
            q                    : (vector) of length n
            P                    : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
            date                 : (string) date for saving the results
            frank_wolfe_variants : (vector) of strings, variants of the Frank-Wolfe algorithm
            off_the_shelves      : (vector) of strings, off-the-shelf methods
        Optional Input:
            dual_comparison      : (boolean) if true, the dual gap is computed
            eps_RDG              : (float) max relative duality gap (deal error)
            eps_RT               : (float) max relative tolerance for off-the-shelf methods
            eps_RE               : (float) max relative error to consider a solution correct
            max_steps            : (integer) stop criterion for Frank Wolfe (max number of steps for Frank Wolfe)
        Output:
            table_results        : (table) results obtained by the methods. They are sorted by the primal error time, 
                and number of steps
            table_solutions      : (table) table of solutions (x_min) for each method
    %}

    numvarargs = length(varargin);
    if numvarargs > 5
        error('myfuns:ComparingMethods:TooManyInputs', ...
            'requires at most 5 optional inputs');
    end

    % Set defaults for optional inputs
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
    Primal_Error = zeros(0,1);
    Converging = zeros(0,1);
    Feasible = zeros(0,1);
    Dual_Error = zeros(0,1);

    Solutions = zeros(n, 0);
    Histories = cell(1, 0);

    % Frank-Wolfe type algorithms
    for i = 1:length(frank_wolfe_variants)
        [x_min, f_min, elapsed_time, num_steps, type, variant, primal_error, dual_error, converging, feasible, history] = ...
            FrankWolfe(Q, q, P, frank_wolfe_variants(i), eps_RDG, eps_RE, max_steps, false, false, date, f_star);
        Method = cat(1, Method, type);
        Variants = cat(1, Variants, variant);
        Minima = cat(1, Minima, f_min);
        Time = cat(1, Time, elapsed_time);
        Steps = cat(1, Steps, num_steps);
        Primal_Error = cat(1, Primal_Error, primal_error);
        Dual_Error = cat(1, Dual_Error, dual_error);
        Converging = cat(1, Converging, converging);
        Feasible = cat(1, Feasible, feasible);
        Solutions = cat(2, Solutions, x_min);
        Histories = cat(2, Histories, history);
    end

    % Quadratic Programming algorithms
    for i = 1:length(off_the_shelves)
        [x_min, f_min, elapsed_time, num_steps, type, variant, primal_error, dual_error, converging, feasible, history] = ...
            QuadraticProgramming(Q, q, P, off_the_shelves(i), max_steps, eps_RT, eps_RE, f_star);
        Method = cat(1, Method, type);
        Variants = cat(1, Variants, variant);
        Minima = cat(1, Minima, f_min);
        Time = cat(1, Time, elapsed_time);
        Steps = cat(1, Steps, num_steps);
        Primal_Error = cat(1, Primal_Error, primal_error);
        Dual_Error = cat(1, Dual_Error, dual_error);
        Converging = cat(1, Converging, converging);
        Feasible = cat(1, Feasible, feasible);
        Solutions = cat(2, Solutions, x_min);
        Histories = cat(2, Histories, history);
    end

    table_results = table(Method, Variants, Minima, Time, Steps, Primal_Error, Dual_Error, Converging, Feasible);
    table_results = sortrows(table_results, {'Primal_Error', 'Time', 'Steps'});

    label_column = string(frank_wolfe_variants);
    label_column = cat(2,label_column,string(off_the_shelves));
    table_solutions = array2table(Solutions, 'VariableNames', label_column);

    colors = ["#0072BD","#D95319","#77AC30","#EDB120","#7E2F8E"];
    iterative_algorithms = cat(2,append('FW-',string(frank_wolfe_variants)),off_the_shelves);
    PlotBenchmark(Histories, f_star, iterative_algorithms, dual_comparison, colors, date)
end