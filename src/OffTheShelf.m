function [x_min, f_min, elapsed_time, num_steps, method,  variant, converging, feasible, duality_gap, history] = OffTheShelf(Q, q, P, algorithm, max_steps, eps_R, varargin)
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

fun = @(x) x'*Q*x + q'*x;
x0 = StartingPoint(P);
Aeq = P;
beq = ones(size(Aeq,1),1);
lb = zeros(size(Aeq,2),1);
ub = inf(size(Aeq,2),1);

history = [];
function stop = outfun(x,optimValues,state)
    stop = false;
    %disp([optimValues.fval,eps_R * max(1,abs(f_star)),optimValues.fval - f_star]);
    %CheckDomain(x,P)
    % Check conditions to see whether optimization should quit
    if isequal(state,'iter')
        history = [history; optimValues.fval];
    end
    if((StoppingCriteria(optimValues.fval, f_star, eps_R) && CheckDomain(x,P)) || optimValues.iteration>=max_steps)
        stop = true;
    end
end

%CheckDomain(x0,P)

void_options = optimoptions(@fmincon,'OutputFcn',@outfun, 'Algorithm', algorithm, ...
    'MaxFunctionEvaluation', inf, ...
    'MaxIterations', inf, ...
    'StepTolerance', 0, ...
    'ConstraintTolerance', 0, ...
    'OptimalityTolerance', 0, ...
    'ObjectiveLimit', -inf, 'Display', 'off');

disp(strcat('Off-the-shelf constrained optimization: ', algorithm))

tic
[x_min, f_min, ~, output] = fmincon(fun, x0, [], [], Aeq, beq, lb, ub, [], void_options);
elapsed_time = toc;

method = algorithm;
num_steps = output.iterations;
%converging = (exitflag == -1);
%feasible = (exitflag ~= -2);
feasible = CheckDomain(x_min,P);
converging = StoppingCriteria(f_min, f_star, eps_R) && feasible;

end