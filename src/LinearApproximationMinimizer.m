function [d, y, duality_gap] = LinearApproximationMinimizer(Q, q, x, indices, partition)
%{
Minimize the linear approximation of the function in the point x; compute the
descent direction and the "residual error"
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
    duality_gap : (float) scalar product between the gradient in x and the descent direction
%}    

% Gradient function
Df = @(x) 2*Q*x + q;

D = Df(x);
n = size(Q, 1);

y = zeros(n, 1);
y(indices) = 1;

for simplex = partition
    % Take the non-zero indices (indices in I_k)
    Ik = simplex{:};
    
    % Compute the argmin of D restricted on Ik
    [~, j] = min(D(Ik));
    
    % Insert 1 at the position j_k
    y(Ik(j)) = 1;
end

% Define the descent direction
d = y - x;

% Scalar product between the descent direction and the gradient in x 
duality_gap = - d' * D;

end