function PlotLineSearch(Q, q, x, d, alpha, alpha_max, i)
    %{
        Plot the tomography
        Input:
            Q         : (matrix) nxn positive semi-definite
            q         : (vector) of length n
            x         : (vector) starting point
            d         : (vector) descent direction
            alpha     : (float) coefficient for the descent direction
            alpha_max : (float) maximum value of alpha
            i         : (integer) iteration number
    %}

    f = @(t) (x+t*d)'*Q*(x+t*d) + q'*(x+t*d);

    t=0;
    y=f(t);
    delta=(2*alpha)/100;
    stop = (2*alpha)*(2*alpha<1) + 1*(2*alpha >=1);
    while(t < stop)
        t=t+delta;
        y=cat(1,y,f(t));
    end

    gcf = figure('Name','Tomography');
    plot(linspace(0,stop,length(y)),y, 'k', 'LineWidth',2)
    title(['From ', 'f(x_{', num2str(i), '})', '  to  ', 'f(x_{', num2str(i+1), '}) with \alpha_{max} = ', num2str(alpha_max)])
    xlabel('alpha')
    ylabel('f(x)')
    hold on
    plot(alpha, f(alpha), '*', 'Color', 'red', 'MarkerSize',12, 'LineWidth',2)
    if alpha_max < stop
        y_lim = ylim;
        plot([alpha_max alpha_max], [y_lim(1) y_lim(2)], '--', 'Color', 'black', 'LineWidth',2)
    end
    hold off
    fontsize(gcf,scale=1.4)
    grid on
end