function [x_min, f_min]= FrankWolfe_visualize(Q, q, P, x_start, eps, eps_ls, line_search, beta)
%{
FrankWolfe computes the minimum of a quadratic function in a constrained convex domain,
    showing a plot for each iteration. 
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    P           : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    x_start     : (vector) start point
    eps         : (float) stop criterion for Frank Wolfe
    eps_ls      : (float) stop criterion for the line search
    line_search : (string) method for line search
    beta        : (float) momentum coefficient
Output:
    x_min       : (vector) argmin of the function
    f_min       : (vector) min of the function
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
    disp('Using momentum')
end

% Force q and x to be column vectors
q = q(:);
x = x_start(:);

% Function f
f = @(x) x'*Q*x + q'*x;
% Gradient function
Df = @(x) 2*Q*x + q;

[K, n] = size(P);

% Start the algorithm
i = 0;
D = Df(x);
y = zeros(n, 1);
fx(1) = f(x);
disp(['it. ', num2str(i), ', f(x) = ', num2str(fx(1))])

figure('Name','Main');
w = waitforbuttonpress;

for k = 1 : K
    % Take the indeces not equal to zero (indeces in I_k)
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
obj = D' * d;
O(1) = obj;

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

% Plot the line search
plotLS(Q, q, x, d, alpha, alphaStart)

% Upgrade the vector x
x = x + alpha * d;

% Evaluate the function at the new point
fx(2) = f(x);

% Save the momentum direction
d_old = y - x;

% Print the results of the first iteration
disp('Gradiente Df(x):')
disp(D)
disp('argmin <y,Df(x)>:')
disp(y)
disp(['<grad, d> = ', num2str(obj)])
disp(['it. ', num2str(i + 1), ' alpha = ', num2str(alpha), ' f(x) = ', num2str(fx(2))])

% New iteration
i = i + 1;

w = waitforbuttonpress;

% Iterate until convergence
while (obj < - eps)
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
    obj = D'*d;
    O(i) = obj;
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
    %   the new point must be inside the triangular with vertices x, d_old,
    %   d_new
    par_momentum = min(beta, 1 - alpha);
    par_momentum = max(0, par_momentum);
    if(beta > 0)
        plotMOMENTUM(Q, q, x, d_old, par_momentum, d, alpha)
    else
        plotLS(Q, q, x, d, alpha, alphaStart)
    end
    momentum = par_momentum * d_old;
    x = x + alpha * d + momentum;
    fx(i+1) = f(x);
    d_old = y - x;
    i = i + 1;
    disp('Gradiente Df(x):')
    disp(D)
    disp('argmin <y,Df(x)>:')
    disp(y)
    disp(['<grad, d> = ', num2str(obj)])
    disp(['it. ', num2str(i), ', alpha = ', num2str(alpha), ', f(x) = ', num2str(fx(i))])
    if (beta > 0)
        disp(['par_momentum = ', num2str(par_momentum)])
    end
    w = waitforbuttonpress;
end

figure('Name','Objective Function');
plot(O, 'ro-')
hold on
plot(fx, 'bo-')
hold off

x_min = x;
f_min = f(x);
