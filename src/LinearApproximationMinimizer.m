function [d, y, duality_gap] = LinearApproximationMinimizer(Q, q, x, indices, partition)
%{
Minimize the linear approximation of the function in the point x; compute the
    descent direction and the "residual error"
INPUT:
    Q         : (matrix) nxn positive semi-definite
    q         : (vector) of length n
    x         : (vettor) point where compute and minimize the linear approximation of the funtion
    indices   : (vector) indices belongind to the simplices with at least two vertices
    partition : (cell-array) each array represents a simplex
OUTPUT:
    d           : (vector) descent direction
    y           : (vector) point that minimizes the dot product with the gradient of the function
        in the point x
    duality_gap : (float) scalar product between the gradient in x and the descent direction
%}    

if nargin < 5
    indices = true(1,size(Q,1));
end

% Gradient function
Df = @(x) 2*Q*x + q;

%{
[~, n] = size(P);
D = Df(x);
DD = kron(D,ones(1,n));
indices = logical(P)';
DD(~indices) = inf;
[~, ind_min] = min(DD, [], 1);
y(ind_min) = 1;
%}
%{
[K, n] = size(P);
y = zeros(n, 1);
D = Df(x);
mask = logical(P)';
for k = 1 : K
    % Extract the gradient components
    [~, j] = min(D(mask(:,k)));
    
    % Insert 1 at the position j_k
    y = y(mask(:,k));
    y(j) = 1;
end
%}

D = Df(x);
n = size(Q, 1);
y = zeros(n, 1);
x(~indices) = 0;

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