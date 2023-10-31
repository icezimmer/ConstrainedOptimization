function [d_a, y_a, duality_gap_a, alpha_max_a] = LinearizationMaximizer(Q, q, x, partition)
    %{
        Maximize the linear approximation of the function in the point x, compute the
        descent direction () and the duality gap
        Input:
            Q             : (matrix) nxn positive semi-definite
            q             : (vector) of length n
            x             : (vector) point where compute and minimize the linear approximation of the funtion
            partition     : (cell-array) each cell contains the indices of the vertices of a simplex
        Output:
            d_a           : (vector) descent direction (opposite of away direction)
            y_a           : (vector) point where the linear approximation is maximized
            duality_gap_a : (scalar) duality gap respect the opposite of away direction
            alpha_max_a   : (scalar) maximum step size respect the opposite of away direction 
    %}    

    % Gradient function
    grad = @(x) 2*Q*x + q;

    grad_x = grad(x);
    n = size(Q, 1);

    y_a = zeros(n, 1);

    % Search the active indices
    active_set = x > 0;

    alpha_max_a = inf;
    for simplex = partition

        % Take the indices in I_k
        I_k = simplex{:};

        % Take the active indices of I_k
        S_k = I_k(active_set(I_k));
        % Compute the argmax of D restricted on Sk
        [~, j_max] = max(grad_x(S_k));     
        % Insert 1 at the position j_max
        y_a(S_k(j_max)) = 1;
        % Avoid division by denominator < 0
        new_alpha_max = x(S_k(j_max)) / max(0, 1-x(S_k(j_max)));
        if new_alpha_max < alpha_max_a
            alpha_max_a = new_alpha_max;
        end

    end

    % Compute the opposite of away direction
    d_a = x - y_a;
    % Compute the duality gap respect the opposite of away direction
    duality_gap_a = - d_a' * grad_x;
end