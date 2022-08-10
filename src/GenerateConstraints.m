function P = GenerateConstraints(n, K, seed)
%{
Construct the matrix P representing the partion of indices {I_k}. Each row represents a simplex.
Input:
    n    : (integer) dimension of the space
    K    : (integer) number of simplices
    seed : (integer) seed for the random generator
Output:
    P    : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
%}

rng(seed)

aux=sort(randperm(n,K));
dim_Simplex=diff(aux);
dim_Simplex(end+1)=n-sum(dim_Simplex);

P = zeros(K, n);
P(1, 1 : dim_Simplex(1)) = ones(1, dim_Simplex(1));
column = dim_Simplex(1);
for row = 2 : K
    P(row, column + 1 : column + dim_Simplex(row)) = ones(1, dim_Simplex(row));
    column = column + dim_Simplex(row);
end

% Randomly shuffle the columns of the matrix P
shuffle = randperm(n);
P = P(:,shuffle);

end