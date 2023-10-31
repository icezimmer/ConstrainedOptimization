function linearly_independent_rows = FindLinearlyIndependentRows(A, k)
    %{
        Find all linearly independent combinations of k rows from A
        Input:
            A : (matrix) The matrix to find linearly independent rows from
            k : (integer) The number of rows to find linearly independent combinations of
        Output:
            linearly_independent_rows : (cell array) A cell array containing all linearly independent combinations 
                of k rows from A
    %}
    % Get the number of rows in A
    n = size(A, 1);
    
    % Compute all combinations of k rows from A
    row_combinations = nchoosek(1:n, k);
    
    % Initialize an empty cell array to store the independent matrices
    linearly_independent_rows = cell(1, size(row_combinations, 1));
    
    % Iterate over each row combination
    for i = 1:size(row_combinations, 1)
        % Extract the rows corresponding to the current combination
        rows = A(row_combinations(i, :), :);
        
        % Check if the rows are linearly independent
        if rank(rows) == k
            % Store the independent matrix in the cell array
            linearly_independent_rows{i} = rows;
        end
    end
    
    % Remove any empty cells from the cell array
    linearly_independent_rows = linearly_independent_rows(~cellfun('isempty', linearly_independent_rows));
end
