function [x_star, f_star] = Oracle(Q, q, P)
    %{
        Computes the constrained minimum of f(x) = x'Qx + q'x with high precision using interior-point-convex
        Input:
            Q      : (matrix) nxn positive semi-definite
            q      : (vector) of length n
            P      : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
        Output:
            x_star : (vector) of length n (solution of task)
            f_star : (float) f(x_star)
    %}
    disp('Pre-compute the optimum to compute the primal error (relative error)')

    H = 2 * Q;
    vec = q;
    Aeq = P;
    beq = ones(size(Aeq,1),1);
    lb = zeros(size(Aeq,2),1);

    options = optimoptions(@quadprog, 'Algorithm', "interior-point-convex", ...
        'MaxIteration', inf, ...
        'OptimalityTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10, ...
        'StepTolerance', 1e-12, ...
        'Display', 'off');

    % Compute the solution using interior-point-convex
    [x_star, f_star] = quadprog(H, vec, [], [], Aeq, beq, lb, [], [], options);
end