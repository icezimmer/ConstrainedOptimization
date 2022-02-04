function [x_min, f_min, elapsed_time, num_steps, converging, duality_gap] = FrankWolfe(Q, q, P, x_start, eps, max_steps, eps_ls, line_search, beta, tomography, curve)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain. 
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    P           : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start     : (vector) start point
    eps         : (float) stop criterion (max duality_gap for Frank Wolfe)
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
    duality_gap  : (float) opposite value given by the scalar product between
        the gradient and the direction 
%}

if isequal(line_search,'NM')
    disp('LineSearch: Newton Method')
elseif isequal(line_search,'LBM')
    disp('LineSearch: Linear Bisection Method')
elseif isequal(line_search,'QBM')
    disp('LineSearch: Quadratic Bisection Method')
elseif isequal(line_search,'Default')
    disp('LineSearch: Default Method')
end

if (beta > 0)
    disp(strcat("Using momentum = ", num2str(beta)))
end

% Force q and x to be column vectors
q = q(:);
x = x_start(:);

% Function f
f = @(x) x'*Q*x + q'*x;

tic
% First iteration
i = 0;

[d, y, gap] = LinearApproximationMinimizer(Q, q, P, x);

% Line search
alpha = LineSearch(Q, q, x, d, eps_ls, i, line_search);

if tomography
    figure('Name','Main');
	w = waitforbuttonpress;
    % Print the f value for the start point
	disp(['it. ', num2str(i), ', f(x) = ', num2str(f(x))])
    % Plot the tomography
    if ~isequal(line_search, 'Default')
        % Plot the line search
        plotLS(Q, q, x, d, alpha, alphaStart, i)
    end
    if isequal(line_search, 'Default')
        plotLS(Q, q, x, d, alpha, 1, i)
    end
    % Print the direction and the alpha step computed at the first iteration
    disp(['<grad, d> = ', num2str(gap), ', alpha = ', num2str(alpha)])
end

if curve
    fx(i+1) = f(x);
    E(i+1) = gap;
end

% Upgrade the vector x
x = x + alpha * d;

% Save the momentum direction
d_old = y - x;

% New iteration
i = i + 1;

% Iterate until convergence
while (gap < - eps && i < max_steps)
    [d, y, gap] = LinearApproximationMinimizer(Q, q, P, x);
    
    alpha = LineSearch(Q, q, x, d, eps_ls, i, line_search);
    
    % Compute the momentum:
    [momentum, momentum_coeff] = Momentum(d_old, alpha, beta);
    
    % Plot tomography
    if tomography
        disp(['it. ', num2str(i), ', f(x) = ', num2str(f(x))])
        w = waitforbuttonpress;
        if ~isequal(line_search, 'Default')
            if(beta > 0)
                plotMOMENTUM(Q, q, x, d_old, momentum_coeff, d, alpha)
            else
                plotLS(Q, q, x, d, alpha, alphaStart, i)
            end
        end
        if isequal(line_search, 'Default')
            if(beta > 0)
                plotMOMENTUM(Q, q, x, d_old, par_momentum, d, alpha)
            else
                plotLS(Q, q, x, d, alpha, 1, i)
            end
        end
        disp(['<grad, d> = ', num2str(gap), ', alpha = ', num2str(alpha)])
    end
    
    % Append new value of the funtion and the duality_gap for the optimization curve
    if curve
        fx(i+1) = f(x);
        E(i+1) = gap;
    end
    
    % Upgrade the point
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
num_steps = i;

% Check the convergence of the algorithm
if(gap < - eps || ~Domain(x_min, P))
    converging = "No";
else
    converging = "Yes";
end

% Final duality_gap (duality gap) 
duality_gap = -E(end);

if tomography
    disp(['it. ', num2str(i), ', f(x) = ', num2str(f_min)])
end

% Plot the optimization curve
if curve
    fx(i+1) = f_min;
    plotCURVE(fx, E, line_search, beta)
end

end