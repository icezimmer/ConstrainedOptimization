function [x_min, f_min, elapsed_time, num_steps, method,  variant, err, converging, feasible, duality_gap, history] = QuadraticProgramming(Q, q, P, algorithm, max_steps, eps_RT, eps_RE, varargin)
%{
Quadratic programming using the built-in function "quadprog" by the optimization-toolbox of MATLAB
Input:
    Q         : (matrix) nxn positive semi-definite
    q         : (vector) of length n
    P         : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    algorithm : (string) name of algorithm for the quadratic programming
    max_steps : (integer) stop criterion (max number of steps for Frank Wolfe)
    f_star    : (float) optimal value of the function in the domain
Output:
    x_min            : (vector) argmin of the function
    f_min            : (vector) min of the function
    elapsed_time     : (float) time elapsed for the computation
    step_size_method : (string)
    num_steps        : (integer) number of steps for the convergence
    converging       : (logical) method converge or not
    feasible         : (logical) solution is feasible or not
    duality_gap      : (float) opposite value of the scalar product between
        the descent direction and the gradient
%}

numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:OffTheShelf:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

if numvarargs == 1
    f_star = varargin{1};
elseif numvarargs < 1
    [~, f_star] = Oracle(Q, q, P);
end

variant = "--";
duality_gap = NaN;

x0 = StartingPoint(P);
H = 2 * Q;
vec = q;
Aeq = P;
beq = ones(size(Aeq,1),1);
lb = zeros(size(Aeq,2),1);

start_options = optimoptions(@quadprog, 'Algorithm', algorithm, ...
    'MaxIteration', max_steps-1, ...
    'OptimalityTolerance', eps_RT, ...
    'ConstraintTolerance', eps_RT, ...
    'StepTolerance', eps_RT, ...
    'Display', 'iter-detailed');

f = @(x) x'*Q*x + q'*x;

disp(strcat('Off-the-shelf constrained optimization: ', algorithm))

tic
[x_min, f_min, ~, output] = quadprog(H, vec, [], [], Aeq, beq, lb, [], x0, start_options);
elapsed_time = toc;
num_steps = output.iterations;
history = zeros(1,num_steps+1);
history(1) = f(x0);
history(end) = f_min;

for step = 0:num_steps-2
    options = optimoptions(start_options, 'MaxIteration', step, 'Display', 'off');
    [~, f_x] = quadprog(H, vec, [], [], Aeq, beq, lb, [], x0, options);
    history(step+2) = f_x;
end

method = algorithm;
feasible = CheckDomain(x_min,P);
[converging, err] = ConvergingError(f_min, f_star, eps_RE);
converging = converging & feasible;

end