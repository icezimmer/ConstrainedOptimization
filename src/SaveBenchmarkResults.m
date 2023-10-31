function SaveBenchmarkResults(table_results, table_solutions, date)
    %{
        Save the results of the benchmark in a .mat file
        Input:
            table_results   : (table) results obtained by the methods. They are sorted by the minimum,
                duality_gap, time and number of steps
            table_solutions : (table) table of solutions (x_min) for all the methods
            date            : (float) date for saving the results 
    %}

    path = fullfile('results', date, 'benchmark_results.mat');

    benchmark_results.table_results = table_results;
    benchmark_results.table_solutions = table_solutions;

    save(path, "benchmark_results");
end