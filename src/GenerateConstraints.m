function P = GenerateConstraints(n, seed)
%{
Construct the matrix P representing the partion of indices {I_k}. Each row represents a simplex.
Input:
    n    : (integer) dimension of the space
    seed : (integer) seed for the random generator
Output:
    P    : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
%}

rng(seed)

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

% Randomly shuffle the columns of the matrix P
shuffle = randperm(n);
P = P(:,shuffle);

end