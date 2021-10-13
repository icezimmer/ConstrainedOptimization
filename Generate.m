function [Q, q, P, xStart] = Generate(n, dim_Ker)
%Funzione che genera Q, q, P 'quasi' random e xStart nel dominio
%INPUT n = dimensione

A = -2 + 4 * rand(n);
A(n - dim_Ker + 1 : end, :) = zeros(dim_Ker, n);
%{
r = zeros(1, 0);
n0 = n;
while (n > 0)
    r0 = randi([1, n]);
    r = [r, r0];
    n = n - r0;
end

r0 = randi([1, n0 - 1]);
A = -2 + 4 * rand(n0);
A(n0 - r0 + 1 : end, :) = zeros(r0, n0);
%}

Q = A' * A;

q = -10 + 20 * rand(n, 1);

q_plus = Q\q;
q_0 = Q * q_plus - q;

r = zeros(1, 0);
n0 = n;
while (n > 0)
    r0 = randi([1, n]);
    r = [r, r0];
    n = n - r0;
end
m = length(r);
P = zeros(m, n);
P(1, 1 : r(1)) = ones(1, r(1));
s = r(1);
for i = 2 : m
    P(i, s + 1 : s + r(i)) = ones(1, r(i));
    s = s + r(i);
end

xStart = zeros(n, 1);
xStart(1 : r(1)) = ones(r(1), 1)/r(1);
s = r(1);
for i = 2 : m
    xStart(s + 1 : s + r(i)) = ones(r(i), 1)/r(i);
    s = s + r(i);
end

disp(['norm(q_0) = ', num2str(norm(q_0))])