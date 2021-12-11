function [Q, q, P, x_start, minima] = Generate(n, dim_ker, min_eig, max_eig, min_q, max_q, zero_q)
%{
Generate randomly the matrix Q, the vector q, the matrix P (representing the partion
    of indices) and the point x_start belonging to the domain.
Input:
    n       : (integer) dimension of the space
    dim_Ker : (integer) dimension of the kernel space
    min_eig : (float) min eigenvalue of the matrix Q
    max_eig : (float) max eigenvalue of the matrix Q
    min_q   : (float) min value of the vector q
    max_q   : (float) max value of the vector q
    zero_q  : (integer) number of zero values in the vector q
Output:
    Q      : (matrix) nxn positive semi-definite with non-zero eigenvalues in the
        range (min_eig, max_eig) 
    q      : (vector) of length n with values in the range (min_q, max_q)
    P      : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    xStart : (vector) start point
    minima : (boolean) 1 iff the function has any minima in R^n space
%}

% Construct the orthogonal matrices U and V
U = orth(randn(n));
V = orth(randn(n));

% Construct the diagonal matrix of the root saquares of
%   the eigenvalues of Q
sv = (abs(sqrt(max_eig) - sqrt(min_eig)) * rand(n - dim_ker, 1)) + min(sqrt(min_eig), sqrt(max_eig));
sv = [sv; zeros(dim_ker,1)];
S = diag(sv);

% Construct the matrix A
A = U * S * V;

% Construct the matrix Q
Q = A' * A;

% Construct the vector q
q = min_q + abs(max_q - min_q) * rand(n - zero_q, 1);
q = [q; zeros(zero_q, 1)];
q = q(randperm(n));

% Ceck if the vector q belong to the image of the matrix Q
minima = (rank(Q) == rank([Q, q]));

% Construct the matrix P representing the partion of indices {I_k}
r = zeros(1, 0);
n0 = n;
while (n0 > 0)
    r0 = randi([1, n0]);
    r = [r, r0];
    n0 = n0 - r0;
end
m = length(r);
P = zeros(m, n);
P(1, 1 : r(1)) = ones(1, r(1));
s = r(1);
for i = 2 : m
    P(i, s + 1 : s + r(i)) = ones(1, r(i));
    s = s + r(i);
end

% Construct the starting point for the Frank Wolfe algorithm
x_start = zeros(n, 1);
x_start(1 : r(1)) = ones(r(1), 1) / r(1);
s = r(1);
for i = 2 : m
    x_start(s + 1 : s + r(i)) = ones(r(i), 1) / r(i);
    s = s + r(i);
end

end