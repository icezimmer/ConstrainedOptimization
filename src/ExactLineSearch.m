function alpha = ExactLineSearch(Q, q, x, d)
%{
Compute the exact line search
Input:
    Q : (matrix) nxn positive semi-definite
    q : (vector) of length n
    x : (vector) start point
    d : (vector) descent direction
Output:
    alpha : (float) step-size value
%}

grad = @(x) 2*Q*x + q;

% if the denominator is zero alpha will be infinity so we'll fix alpha to 1
alpha = - ((d'*grad(x)) / (2*d'* Q *d));

if alpha > 1
    alpha = 1;
end


end