function alpha = LineSearchNM(Q, q, x, d, alpha_start, eps)
%{
Compute the line search by the Newton method
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    x           : (vector) start point
    d           : (vector) descent direction
    alpha_start : (float) start value for alpha
    eps         : (float) stop criterion
Output:
    alpha       : (float) step-size value
%}

dPhi = @(alpha) (2 * Q * (x + alpha * d) + q)' * d;
ddPhi = 2 * d' * Q * d;

alpha = alpha_start;

% Apply the Newton method only if the denominator is non-zero
if(norm(ddPhi) > eps)
    while(norm(dPhi(alpha)) > eps && alpha <= 1 + eps)
        alpha = alpha - dPhi(alpha)/ddPhi;
    end
end

if(alpha > 1)
    alpha = 1;
end
if(alpha < 0)
    alpha = 0;
end

end