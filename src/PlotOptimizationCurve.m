function PlotOptimizationCurve(fx, E, line_search)
%{
Plot the optimization curve
Input:
    fx          : (vector) function values sequence
    E           : (vector) error values sequence
    line_search : (string) line search method
%}

figure('Name', strcat(line_search));

tiledlayout(2,1)

nexttile
plot(fx, 'bo-')
title('Primal optimization')
xlabel('step')
ylabel('f(x)')

nexttile
plot(E, 'ro-')
title('Dual optimization')
xlabel('step')
ylabel('duality gap')

end