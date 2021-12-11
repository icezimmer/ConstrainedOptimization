function alpha = LineSearchLBM(Q, q, x, d, alpha_start, eps)
%{
Compute the line search by the linear bisection method
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

alphaMinus = 0;
alpha = alpha_start;
alphaPlus = alpha_start;

% Until the gradient value and the range are too small
while(norm(dPhi(alpha)) > eps && norm(alphaPlus - alphaMinus) > eps)
    alpha = (alphaPlus + alphaMinus) / 2;
    if(dPhi(alpha) < 0)
        alphaMinus = alpha;
    else
        alphaPlus = alpha;
    end
end

end