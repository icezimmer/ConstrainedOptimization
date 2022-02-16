function [x_min, f_min, elapsed_time, method, step_size_method, num_steps, converging, feasible, duality_gap] = DirectComputation(Q, q, P)
%{
Input:
    Q : (matrix) nxn positive semi-definite
    q : (vector) of length n
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

% Function f
f = @(x) x'*Q*x + q'*x;

[n, ~] = size(Q);

x_min = NaN(n, 1);
f_min = NaN;
elapsed_time = NaN;
method = "Direct";
step_size_method = "--";
num_steps = NaN;
duality_gap = NaN;

converging = false;
feasible = false;

if (rank(Q) == rank([Q, q]))
    tic
    x_min = (-2 * Q) \ q;
    f_min = f(x_min);
    converging = true;
    feasible = Domain(x_min, P);
    elapsed_time = toc;
end

end

