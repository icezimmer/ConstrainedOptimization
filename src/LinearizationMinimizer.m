function [d, y, duality_gap] = LinearizationMinimizer(Q, q, x, partition)
    %{
        Minimize the linear approximation of the function in the point x, compute the
        descent direction and the duality gap
        Input:
            Q         : (matrix) nxn positive semi-definite
            q         : (vector) of length n
            x         : (vector) point where compute and minimize the linear approximation of the funtion
            partition : (cell-array) each cell contains the indices of the vertices of a simplex
        Output:
            d           : (vector) descent direction
            y           : (vector) point where the linear approximation is minimized
            duality_gap : (float) opposite of scalar product between the gradient in x and the descent direction
    %}    

    % Gradient function
    grad = @(x) 2*Q*x + q;

    grad_x = grad(x);
    n = size(Q, 1);

    y = zeros(n, 1);

    for simplex = partition
        % Take the indices in I_k
        I_k = simplex{:};
        
        % Compute the argmin of D restricted on I_k
        [~, j] = min(grad_x(I_k));
        
        % Insert 1 at the position j_k
        y(I_k(j)) = 1;
    end

    % Compute the descent direction
    d = y - x;

    % Compute the duality gap
    duality_gap = - d' * grad_x;
end