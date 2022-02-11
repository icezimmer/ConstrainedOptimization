function alpha_start = StartLineSearch(Q, q, x, d, eps)
%{
Initialize the parameter alpha for the line search
Input:
    Q   : (matrix) nxn positive semi-definite
    q   : (vector) of length n
    x   : (vector) start point
    d   : (vector) descent direction
    eps : (float) stop criterion
Output:
    alphaStart : (float) start value for alpha
%}

dPhi = @(alpha) (2 * Q * (x + alpha * d) + q)' * d;

alpha_start = 0.1;

% Until the gradient is negative and alpha is grater than one
while(dPhi(alpha_start) <= -eps && alpha_start <= 1 + eps)
    alpha_start = 2*alpha_start;
end

if (alpha_start > 1)
    alpha_start = 1;
end

end