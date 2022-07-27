function [x_min, f_min, elapsed_time, method, step_size_method, num_steps, converging, feasible, duality_gap] = QuadraticProgramming(Q, q, P, algorithm, max_steps)
%{
Quadratic programming using the built-in function "quadprog" by the optimization-toolbox of MATLAB
Input:
    Q         : (matrix) nxn positive semi-definite
    q         : (vector) of length n
    P         : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    algorithm : (string) name of algorithm for the quadratic programming
    max_steps : (integer) stop criterion (max number of steps for Frank Wolfe)
Output:
    x_min        : (vector) argmin of the function
    f_min        : (vector) min of the function
    elapsed_time : (float) time elapsed for the computation
    num_steps    : (integer) number of steps for the convergence
    converging   : (logical) method converge or not
    feasible     : (logical) solution is feasible or not
    duality_gap  : (float) opposite value of the scalar product between
        the descent direction and the gradient
%}

disp(strcat('Quadratic Programming by the optimization toolbox: ', algorithm))

if nargin < 4 % no algorithm, max_steps
    algorithm = "interior-point-convex";
    max_steps = 1000;
elseif nargin == 4 % no max_steps
    max_steps = 1000;
end

step_size_method = "--";
duality_gap = NaN;

H = 2 * Q;
f = q;
[K, n] = size(P);
Aeq = P;
beq = ones(K,1);
lb = zeros(n,1);

x0 = StartingPoint(P);

options = optimoptions(@quadprog, 'Algorithm', algorithm, 'MaxIterations', max_steps, 'Display', 'off');

tic
[x_min, f_min, exitflag, output] = quadprog(H, f, [], [], Aeq, beq, lb, [], x0, options);
elapsed_time = toc;

method = algorithm;
num_steps = output.iterations;
converging = (exitflag > 0);
feasible = (exitflag ~= -2);

end

