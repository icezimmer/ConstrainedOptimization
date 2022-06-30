function SaveResults(table_results, table_solutions, date)
%{
Print the results
Input:
    table_results   : (table) results obtained by the methods. They are sorted by the minimum,
        duality_gap, time and number of steps
    table_solutions : (table) table of solutions (x_min) for all the methods
    date            : (float) date for saving the results 
%}

path1 = fullfile('results', date, 'comparison.mat');
path2 = fullfile('results', date, 'solutions.mat');
save(path1, "table_results");
save(path2, "table_solutions");

end