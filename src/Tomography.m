function Tomography(Q, q, x, d, alpha, alpha_max, i, duality_gap)
    %{
        Show info and plot tomography for each step of the Frank-Wolfe
        Input:
            Q           : (matrix) nxn positive semi-definite
            q           : (vector) of length n
            x           : (vector) starting point
            d           : (vector) descent direction
            alpha       : (float) step size
            alpha_max   : (float) maximum step size
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
    PlotLineSearch(Q, q, x, d, alpha, alpha_max, i)


    % Print the direction and the alpha step computed at the first iteration
    disp(['Duality Gap = ', num2str(duality_gap), ', alpha = ', num2str(alpha), ', alpha_max = ', num2str(alpha_max)])
end