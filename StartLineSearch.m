function alphaStart = StartLineSearch(Q, q, x, d, eps)
dPhi = @(alpha) (2 * Q * (x + alpha * d) + q)' * d;

alphaStart = 0.1;
while(dPhi(alphaStart) <= -eps && alphaStart <= 1 + eps) %fintanto che si ha decrescita
    alphaStart = 2*alphaStart;
end

if (alphaStart > 1)
    alphaStart = 1;
end
