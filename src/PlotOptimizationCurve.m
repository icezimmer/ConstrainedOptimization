function PlotOptimizationCurve(fx, f_star, E, line_search, date)
%{
Plot the optimization curve
Input:
    fx          : (vector) function values sequence
    f_star      : (float) optimal value of the function in the domain
    E           : (vector) error values sequence
    line_search : (string) line search method
    date        : (string) date for saving figures
%}

gap_R = (fx - f_star) / abs(f_star);

gcf = figure('Name', strcat('curve_FW_', line_search));
relative_gap = semilogy(gap_R, 'ko-','DisplayName','Relative Gap');
hold on
duality_gap = semilogy(E, 'ro-','DisplayName','Duality Gap');
hold off

title('Error plot')
xlabel('step')
xlabel('error')
legend([relative_gap, duality_gap], 'Location','best')

saveas(gcf, fullfile('results', date, strcat('curve_FW_', line_search, '.png')))

end