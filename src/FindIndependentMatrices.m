function independent_matrices = FindIndependentMatrices(A, k)
    % Get the number of rows in A
    n = size(A, 1);
    
    % Compute all combinations of k rows from A
    row_combinations = nchoosek(1:n, k);
    
    % Initialize an empty cell array to store the independent matrices
    independent_matrices = cell(1, size(row_combinations, 1));
    
    % Iterate over each row combination
    for i = 1:size(row_combinations, 1)
        % Extract the rows corresponding to the current combination
        rows = A(row_combinations(i, :), :);
        
        % Check if the rows are linearly independent
        if rank(rows) == k
            % Store the independent matrix in the cell array
            independent_matrices{i} = rows;
        end
    end
    
    % Remove any empty cells from the cell array
    independent_matrices = independent_matrices(~cellfun('isempty', independent_matrices));
end
