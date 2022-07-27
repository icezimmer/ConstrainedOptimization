function Tomography(Q, q, x, d, alpha, i, duality_gap)
%{
Give info about the computation and plot the tomography in agiven iteration
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    x           : (vector) starting point
    d           : (vector) descent direction
    alpha       : (float) coefficient for the descent direction
    i           : (integer) iteration number
    duality_gap : (float) opposite value of the scalar product between
        the descent direction and the gradient
%}

% Function f
f = @(x) x'*Q*x + q'*x;

if i == 0
    disp('Press a keyboard button for each iteration')
end

waitforbuttonpress;

% Print the f value for the start point
disp(['it. ', num2str(i), ', f(x) = ', num2str(f(x))])

% Plot the line search
PlotLineSearch(Q, q, x, d, alpha, i)


% Print the direction and the alpha step computed at the first iteration
disp(['Duality Gap = ', num2str(duality_gap), ', alpha = ', num2str(alpha)])

end