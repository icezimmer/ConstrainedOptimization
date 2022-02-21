function [Q, q, P, x_start, minima, n, dim_ker, spectral_radius, K_plus, date] = GenerateFromExternalData(matrix_file, seed)
%{
Load the matrix Q, and generate the vector q, the matrix P (representing the partion
    of indices) and the point x_start belonging to the domain.
Input:
    matrix_file : (string) file .mtx, where the matrix is loaded
    seed        : (integer) seed for the random generator
Output:
    Q      : (matrix) nxn positive semi-definite with non-zero eigenvalues in the
        range (min_eig, max_eig) 
    q      : (vector) of length n with values in the range (min_q, max_q)
    P      : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    xStart : (vector) starting point
    minima          : (logical) true iff the function has any minima in R^n space
    spectral_radius : (float) spectral radius of the matrix Q
    K_plus          : (integer) number of simplices with at least 2 vertices
    date            : (float) date for saving parameters, figures and results 
%}

disp('Generating the instance')

% Initialize the random seed
rng(seed)

% Load the matrix Q from the external_data directory
path = fullfile('external_data', matrix_file);
Q = mmread(path);

% Compute the spectral radius of Q
spectral_radius = eigs(Q, 1, 'lm');

% Construct the sparse vector q (same density of the matrix Q)
[n, ~] = size(Q);
density = nnz(Q)/numel(Q);
q = sprandn(n, 1, density);

% Ceck if the vector q belong to the image of the matrix Q
minima = (rank(full(Q)) == rank(full([Q, q])));

% Ceck if the vector q belong to the image of the matrix Q
dim_ker = n - rank(full(Q));

% Construct the matrix P representing the partion of indices {I_k}
P = GenerateConstraints(n, seed);

% Compute the number of simplices with at least 2 vertices
K_plus = sum(sum(P,2) >1);

% Construct the starting point for the Frank-Wolfe algorithm
x_start = zeros(n, 1);
[K, ~] = size(P);
for k = 1 : K
    i = find(P(k,:), 1);
    x_start(i) = 1;
end

formatOut = 'mm_dd_yy_HH_MM_SS';
date = datestr(now,formatOut);
mkdir(fullfile('results',date));

end