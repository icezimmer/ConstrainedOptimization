function [Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, varargin)
%{
Generate randomly the matrix Q, the vector q, the matrix P (representing the partion
of indices) and the point x_start belonging to the domain.
Input:
    n               : (integer) dimension of the space
    K               : (integer) number of simplices
    dim_Ker         : (integer) dimension of the kernel space
    spectral_radius : (float) max eigenvalue of the matrix Q
    density         : (float) density of the matrix Q
    min_q           : (float) min value of the vector q
    max_q           : (float) max value of the vector q
    zero_q          : (integer) number of zero values in the vector q
    seed            : (integer) seed for the random generator
Output:
    Q      : (matrix) nxn positive semi-definite with non-zero eigenvalues in the
        range (min_eig, max_eig) 
    q      : (vector) of length n with values in the range (min_q, max_q)
    P      : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    K_plus : (integer) number of simplices with at least 2 vertices
    K_avg  : (integer) average of dimensions of the simplices with at least 2 vertices
    date   : (float) date for saving parameters, figures and results 
%}

numvarargs = length(varargin);
if numvarargs > 7
    error('myfuns:GenerateInstance:TooManyInputs', ...
        'requires at most 7 optional inputs');
end

% set defaults for optional inputs
optargs = {randi([1,n]),false,0,10,1,0,0};
optargs(1:numvarargs) = varargin;
[K, force_non_point_simplices, dim_Ker, spectral_radius, density, actv, seed] = optargs{:};

if K<1 || K > n
    error("Number of simplices K must be an intenger >= 1 and <= n")
end

if dim_Ker<0 || dim_Ker>n
    error("Dimension of Ker must be an intenger >= 0 and <= n")
end

if density<0 || density>1
    error("Density must be >= 0 and <= 1")
end

if actv<0 || actv>1
    error("The fracion of active constraints must be in [0,1]")
end

disp('Generating the instance')

% Initialize the random seed
rng(seed)

% Vector of eigenvalues of Q
if dim_Ker == n
    rc = zeros(1, n);
else
    rc = [abs(spectral_radius) * [1, rand(1, n-dim_Ker-1)], zeros(1, dim_Ker)];
end

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

% Construct the matrix P representing the partition of indices {I_k}
P = GenerateConstraints(n, K, force_non_point_simplices, seed);

% Construct the vector q
% if norm_q > 0
%     q = 0.5 - rand(n, 1);
%     q = norm_q * q/norm(q);
% elseif norm_q == 0
%     q = zeros(n,1);
% else
%     error("Norm of vector q must be >= 0")
% end
[~, q] = ForceSolution(Q, P, actv);

% Compute the number of simplices with at least 2 vertices
K_plus = sum(sum(P,2) >1);

% Compute the average dimension of simplices among simplices with at least 2 vertices
if K_plus == 0
    K_avg = 1;
else
    K_avg = (n-(K-K_plus))/K_plus;
end

date = string(datetime('now','TimeZone','local','Format','d-MMM-y_HH:mm:ss'));
if ~exist(fullfile('results',date), 'dir')
    mkdir(fullfile('results',date));
end

end