function alpha = LineSearchQBM(Q, q, x, d, alphaStart, eps)
dPhi = @(alpha) (2 * Q * (x + alpha * d) + q)' * d;

alphaMinus = 0;
alpha = alphaStart;
alphaPlus = alphaStart;

%se l'intervallo di ricerca è troppo piccolo termina
while(norm(dPhi(alpha)) > eps && norm(alphaPlus - alphaMinus) > eps)
    alpha = (alphaMinus*dPhi(alphaPlus) - alphaPlus*dPhi(alphaMinus)) / (dPhi(alphaPlus) - dPhi(alphaMinus));
    if(dPhi(alpha) < 0)
        alphaMinus = alpha;
    else
        alphaPlus = alpha;
    end
end