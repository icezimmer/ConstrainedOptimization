function [x_min, f_min, elapsed_time, num_steps, converging, error] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, line_search, beta, tomography, curve)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain. 
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    P           : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start     : (vector) start point
    eps         : (float) stop criterion (max error for Frank Wolfe)
    max_steps   : (integer) stop criterion (max number of steps for Frank Wolfe)
    eps_ls      : (float) stop criterion for the line search
    line_search : (string) method for line search
    beta        : (float) momentum coefficient
    tomography  : (logical) plot or not the tomography for each step
    curve       : (logical) plot or not the optimization curve
Output:
    x_min        : (vector) argmin of the function
    f_min        : (vector) min of the function
    elapsed_time : (float) time elapsed for the computation
    num_steps    : (integer) number of steps for the convergence
    converging   : (string) method converge or not
    error        : (float) absolute value given by the scalar product between
        the gradient and the direction 
%}

if isequal(line_search,'NM')
    disp('LineSearch Newton Method')
elseif isequal(line_search,'LBM')
    disp('LineSearch Linear Bisection Method')
elseif isequal(line_search,'QBM')
    disp('LineSearch Quadratic Bisection Method')
else
    disp('LineSearch Trivial Method')
end

if (beta > 0)
    disp(strcat("Using momentum = ", num2str(beta)))
end

% Force q and x to be column vectors
q = q(:);
x = x_start(:);

% Function f
f = @(x) x'*Q*x + q'*x;
% Gradient function
Df = @(x) 2*Q*x + q;

[K, n] = size(P);

tic
% First iteration
i = 1;
fx(i) = f(x);

D = Df(x);
y = zeros(n, 1);

for k = 1 : K
    % Take the non-zero indices (indeces in I_k)
    Ik = find(P(k,:) == 1);
    % Extract the gradient components
    Dk = D(Ik);
    % Compute the argmin of D_k
    [~, j] = min(Dk);
    % Take the j-th idex of I_k
    jk = Ik(j);
    % Insert 1 at the position j_k
    y(jk) = 1;
end
% Define the descent direction
d = y - x;

% Scalar product between the gradient in x and the descent direction
object = D' * d;
E(i) = object;

if tomography
	% Print the first iteration
	disp(['it. ', num2str(i), ', f(x) = ', num2str(fx(i))])
	disp(['<grad, d> = ', num2str(object)])
	figure('Name','Main');
	w = waitforbuttonpress;
end

% Line search
if isequal(line_search,'LBM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchLBM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end  
elseif isequal(line_search,'QBM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchQBM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end 
elseif isequal(line_search, 'NM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchNM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end  
else
    alpha = 2/(i + 2);
end

if tomography
	% Plot the line search
	plotLS(Q, q, x, d, alpha, alphaStart)
end

% Upgrade the vector x
x = x + alpha * d;

% Save the momentum direction
d_old = y - x;

% New iteration
i = i + 1;

% Iterate until convergence
while (object < - eps && i <= max_steps)
	fx(i) = f(x);
    D = Df(x);
    y = zeros(n, 1);
    for k = 1 : K
        Ik = find(P(k,:) == 1);
        Dk = D(Ik);
        [~, j] = min(Dk);
        jk = Ik(j);
        y(jk) = 1;
    end
    d = y - x;
    object = D'*d;
    E(i) = object;
	if tomography
        disp(['it. ', num2str(i), ', alpha = ', num2str(alpha), ', f(x) = ', num2str(fx(i))])
    	disp(['<grad, d> = ', num2str(object)])
    	if (beta > 0)
        	disp(['par_momentum = ', num2str(par_momentum)])
    	end
    	w = waitforbuttonpress;
	end
    if isequal(line_search,'LBM')
        alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
        if (alphaStart <= 1)
            alpha = LineSearchLBM(Q, q, x, d, alphaStart, eps_ls);
        else
            alpha = 1;
        end  
    elseif isequal(line_search,'QBM')
        alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
        if (alphaStart <= 1)
            alpha = LineSearchQBM(Q, q, x, d, alphaStart, eps_ls);
        else
            alpha = 1;
        end 
    elseif isequal(line_search, 'NM')
        alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
        if (alphaStart <= 1)
            alpha = LineSearchNM(Q, q, x, d, alphaStart, eps_ls);
        else
            alpha = 1;
        end  
    else
        alpha = 2/(i + 2);
    end
    % Compute the momentum:
    %   the new point must be inside the triangular with vertices x, d_old 
    %   and d_new
    par_momentum = min(beta, 1 - alpha);
    par_momentum = max(0, par_momentum);
	if tomography
		if(beta > 0)
            plotMOMENTUM(Q, q, x, d_old, par_momentum, d, alpha)
    	else
        	plotLS(Q, q, x, d, alpha, alphaStart)
    	end
	end
    momentum = par_momentum * d_old;
    x = x + alpha * d + momentum;
    d_old = y - x;
    i = i + 1;
end

% Elapsed time for computation
elapsed_time = toc;
% Minimum belonging to the domain
x_min = x;
% Minimum value of the function
f_min = f(x);
% Number of steps
num_steps = i - 1;

% Control the convergence of the algorithm
if(num_steps >= max_steps && object < - eps && Domain(x_min, P))
    converging = "No";
else
    converging = "Yes";
end

% Final error 
error = abs(E(end));

% Plot the optimizaion curve
if curve
    figure('Name', strcat(line_search, " with momentum = ", num2str(beta)));
    plot(E, 'ro-')
    hold on
    fx(i) = f_min;
    plot(fx, 'bo-')
    hold off
end

end