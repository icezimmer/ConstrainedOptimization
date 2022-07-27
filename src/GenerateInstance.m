function [Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, K, dim_Ker, spectral_radius, density, min_q, max_q, zero_q, seed)
%{
Generate randomly the matrix Q, the vector q, the matrix P (representing the partion
    of indices) and the point x_start belonging to the domain.
Input:
    n               : (integer) dimension of the space
    seed            : (integer) seed for the random generator
    dim_Ker         : (integer) dimension of the kernel space
    spectral_radius : (float) max eigenvalue of the matrix Q
    density         : (float) density of the matrix Q
    min_q           : (float) min value of the vector q
    max_q           : (float) max value of the vector q
    zero_q          : (integer) number of zero values in the vector q
Output:
    Q               : (matrix) nxn positive semi-definite with non-zero eigenvalues in the
        range (min_eig, max_eig) 
    q               : (vector) of length n with values in the range (min_q, max_q)
    P               : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    xStart          : (vector) starting point
    minima          : (logical) true iff the function has any minima in R^n space
    spectral_radius : (float) spectral radius of the matrix Q
    K_plus          : (integer) number of simplices with at least 2 vertices
    date            : (float) date for saving parameters, figures and results 
%}

disp('Generating the instance')

if nargin < 2
    K = ceil(0.01*n);
    dim_Ker = 0;
    spectral_radius = 10;
    density = 1;
    min_q = -5;
    max_q = 5;
    zero_q = 0;
    seed = 0;
elseif nargin == 2
    dim_Ker = 0;
    spectral_radius = 10;
    density = 1;
    min_q = -5;
    max_q = 5;
    zero_q = 0;
    seed = 0;
elseif nargin == 3
    spectral_radius = 10;
    density = 1;
    min_q = -5;
    max_q = 5;
    zero_q = 0;
    seed = 0;
elseif nargin == 4
    density = 1;
    min_q = -5;
    max_q = 5;
    zero_q = 0;
    seed = 0;
elseif nargin == 5
    min_q = -5;
    max_q = 5;
    zero_q = 0;
    seed = 0;
elseif nargin == 6
    max_q = 5;
    zero_q = 0;
    seed = 0;
elseif nargin == 7
    zero_q = 0;
    seed = 0;
elseif nargin == 8
    seed = 0;
end

if K<1 || K > n
    error("Number of simplices K must be an intenger >= 1 and <= n")
end

% Initialize the random seed
rng(seed)

% Vector of eigenvalues of Q
rc = [abs(spectral_radius) * [1, rand(1, n-dim_Ker-1)], zeros(1, dim_Ker)];

if density == 1
    sv = sqrt(rc);
    S = diag(sv);
    U = orth(rand(n));
    V = orth(rand(n));
    A = U * S * V;
    Q = A' * A;
else
    Q = sprandsym(n, density, rc);
end

% Construct the vector q
q = min_q + abs(max_q - min_q) * rand(n - zero_q, 1);
q = [q; zeros(zero_q, 1)];
q = q(randperm(n));

% Construct the matrix P representing the partition of indices {I_k}
P = GenerateConstraints(n, K, seed);

% Compute the number of simplices with at least 2 vertices
K_plus = sum(sum(P,2) >1);

% Compute the average dimension of simplices among simplices with at least 2 vertices
if K_plus == 0
    K_avg = 1;
else
    K_avg = (n-(K-K_plus))/K_plus;

formatOut = 'mm_dd_yy_HH_MM_SS';
date = datestr(now,formatOut);
mkdir(fullfile('results',date));

end