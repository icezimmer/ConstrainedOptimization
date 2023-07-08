function [x_star, f_star] = Oracle(Q, q, P)

disp('Pre-compute the optimum for the stop condition')

% Compute the minimum with off-the-shelf method
H = 2 * Q;
vec = q;
Aeq = P;
beq = ones(size(Aeq,1),1);
lb = zeros(size(Aeq,2),1);

options = optimoptions(@quadprog, 'Algorithm', "interior-point-convex", ...
    'MaxIteration', inf, ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', eps, ...
    'StepTolerance', 1e-10, ...
    'Display', 'off');

[x_star, f_star] = quadprog(H, vec, [], [], Aeq, beq, lb, [], [], options);

end