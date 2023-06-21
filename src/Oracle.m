function f_star = Oracle(Q, q, P)

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
    'StepTolerance', eps, ...
    'Display', 'off');

[~, f_star] = quadprog(H, vec, [], [], Aeq, beq, lb, [], [], options);

%{
Aeq = P;
beq = ones(size(Aeq,1),1);
lb = zeros(size(Aeq,2),1);
ub = inf(size(Aeq,2),1);


best_options = optimoptions(@fmincon, 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluation', inf, ...
    'MaxIterations', inf, ...
    'StepTolerance', 0, ...
    'ConstraintTolerance', 0, ...
    'OptimalityTolerance',1e-6, ...
    'ObjectiveLimit', -inf, 'Display', 'off');

[~, f_star] = fmincon(fun, x0, [], [], Aeq, beq, lb, ub);
%}

end