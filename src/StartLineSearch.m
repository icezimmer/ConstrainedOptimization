function alphaStart = StartLineSearch(Q, q, x, d, eps)
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

alphaStart = 0.1;

% Until the gradient is negative and alpha is grater than one
while(dPhi(alphaStart) <= -eps && alphaStart <= 1 + eps)
    alphaStart = 2*alphaStart;
end

if (alphaStart > 1)
    alphaStart = 1;
end

end