function PrintResults(table_results, table_solutions, date)
%{
Print the results
Input:
    table_results   : (table) results obtained by the methods. They are sorted by the minimum,
        duality_gap, time and number of steps
    table_solutions : (table) table of solutions (x_min) for all the methods
    date            : (float) date for saving the results 
%}

writetable(table_results, fullfile('results', date, 'comparison.txt'), 'Delimiter', 'tab')

writetable(table_solutions, fullfile('results', date, 'solutions.txt'), 'Delimiter', 'tab')

end