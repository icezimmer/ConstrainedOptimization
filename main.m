function [x_min, solution, grid_search] = main()
%{
Search the minimum of a funtion f = x'*Q*x + q'*x
%}

% Space dimension
n = 10;
% Kernel dimension of the matrix Q
dim_ker = 3;
% Minimum eigenvalue of the matrix Q
min_eig = 2;
% Maximum eigenvalue of the matrix Q
max_eig = 10;
% Minimum value of the vector q
min_q = 2;
% Maximum value of the vector q
max_q = 3;
% Number of zero in the vector q
zero_q = 3;

% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, x_start, minima] = Generate(n, dim_ker, min_eig, max_eig, min_q, max_q, zero_q);

% Stop criterion for the Frank-Wolfe method
eps = 0.1;
% Stop criterion for the line search method
eps_ls = 0.01;
% Minimum value of momentum coefficient
left = 0;
% Maximum value of momentum coefficient
right = 0;
% Width of the sub-intervals for the momentum coefficient
delta = 0;

% Search the best solution
[x_min, solution, grid_search] = GridSearch(Q, q, P, x_start, eps, eps_ls, left, right, delta);