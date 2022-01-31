function [d, y, gap] = LinearApproximationMinimizer(Q, q, P, x)
%{
Minimize the linear approximation of the function in the point x; compute the
    descent direction and the "residual error"
INPUT:
    Q : (matrix) nxn positive semi-definite
    q : (vector) of length n
    P : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x : (vettor) point where compute and minimize the linear approximation of the funtion
OUTPUT:
    d   : (vector) descent direction
    y   : (vector) point that minimizes the dot product with the gradient of the function
        in the point x
    gap : (float) scalar product between the gradient in x and the descent direction
%}

% Gradient function
Df = @(x) 2*Q*x + q;

D = Df(x);

[K, n] = size(P);

y = zeros(n, 1);

for k = 1 : K
    % Take the non-zero indices (indeces in I_k)
    Ik = find(P(k,:) == 1);
    
    % Extract the gradient components
    Dk = D(Ik);
    
    % Compute the argmin of D_k
    [~, j] = min(Dk);
    
    % Take the j-th idex of I_k
    jk = Ik(j);
    
    % Insert 1 at the position j_k
    y(jk) = 1;
end

% Define the descent direction
d = y - x;

% Scalar product between the gradient in x and the descent direction
gap = D' * d;

end