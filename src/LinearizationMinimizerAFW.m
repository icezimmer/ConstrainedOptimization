function [d, y, duality_gap, d_a, y_a, duality_gap_a, alpha_max] = LinearizationMinimizerAFW(Q, q, x, indices, partition)
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
y(indices) = 1;

y_a = zeros(n, 1);
% Search the active indices
active_set = x > 0;
y_a(indices & active_set) = 1;

alpha_max = inf;
for simplex = partition

    % Take the indices in I_k
    I_k = simplex{:};
    % Compute the argmin of D restricted on Ik
    [~, j_min] = min(grad_x(I_k));
    % Insert 1 at the position j_min
    y(I_k(j_min)) = 1;

    % Take the active indices of I_k
    S_k = I_k(active_set(I_k));
    % Compute the argmax of D restricted on Sk
    [~, j_max] = max(grad_x(S_k));     
    % Insert 1 at the position j_max
    y_a(S_k(j_max)) = 1;
    % Avoid division by denominator < 0
    new_alpha_max = x(S_k(j_max)) / max(0, 1-x(S_k(j_max)));
    if new_alpha_max < alpha_max
        alpha_max = new_alpha_max;
    end

end

% Compute the descent direction
d = y - x;
% Compute the duality gap respect the descent direction
duality_gap = - d' * grad_x;

% Compute the opposite of away direction
d_a = x - y_a;
% Compute the duality gap respect the opposite of away direction
duality_gap_a = - d_a' * grad_x;

end