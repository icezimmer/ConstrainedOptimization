function alpha = LineSearchNM(Q, q, x, d, alpha, eps)
dPhi = @(alpha) (2 * Q * (x + alpha * d) + q)' * d;
ddPhi = 2 * d' * Q * d;

%posso applicare il NM solo se il denominatore e' diverso da zero
if(norm(ddPhi) > eps)
    while(norm(dPhi(alpha)) > eps && alpha <= 1 + eps)
        alpha = alpha - dPhi(alpha)/ddPhi;
    end
else
    alpha = LineSearchLBM(Q, q, x, d, alpha, eps);
end

if(alpha > 1)
    alpha = 1;
end
if(alpha < 0)
    alpha = 0;
end