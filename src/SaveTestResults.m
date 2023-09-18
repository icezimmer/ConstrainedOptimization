function SaveTestResults(x_min, f_min, elapsed_time, num_steps, method, variant, primal_error, dual_error, converging, feasible, history, date)
%{
Print the results
Input:
    table_results   : (table) results obtained by the methods. They are sorted by the minimum,
        duality_gap, time and number of steps
    table_solutions : (table) table of solutions (x_min) for all the methods
    date            : (float) date for saving the results 
%}

path = fullfile('results', date, 'test_results.mat');

test_results.x_min = x_min;
test_results.f_min = f_min;
test_results.elapsed_time = elapsed_time;
test_results.num_steps = num_steps;
test_results.method = method;
test_results.variant = variant;
test_results.primal_error = primal_error;
test_results.dual_error = dual_error;
test_results.converging = converging;
test_results.feasible = feasible;
test_results.history = history;

save(path, "test_results");

end