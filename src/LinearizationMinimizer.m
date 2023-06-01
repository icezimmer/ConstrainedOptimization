function [d, y, duality_gap] = LinearizationMinimizer(Q, q, x, partition)
%{
Minimize the linear approximation of the function in the point x, compute the
descent direction and the duality gap
INPUT:
    Q         : (matrix) nxn positive semi-definite
    q         : (vector) of length n
    x         : (vector) point where compute and minimize the linear approximation of the funtion
    indices   : (vector) indices belonging to the simplices with only one vertex
    partition : (cell-array) each array represents a simplex with at least two vertices
OUTPUT:
    d           : (vector) descent direction
    y           : (vector) point that minimizes the dot product with the gradient of the function
        in the point x
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