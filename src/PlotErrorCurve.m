function PlotErrorCurve(fx, f_star, E, line_search, date)
%{
Plot the error curve
Input:
    fx          : (vector) function values sequence
    f_star      : (float) optimal value of the function in the domain
    E           : (vector) error values sequence
    line_search : (string) line search method
    date        : (string) date for saving figures
%}

gap_A = (fx - f_star);
gap_R = gap_A / abs(f_star);

gcf = figure('Name', strcat('curve_FW_', line_search));
absolute_gap = semilogy(gap_A, 'bo-','DisplayName','Absolute Error');
hold on
relative_gap = semilogy(gap_R, 'ro-','DisplayName','Relative Error');
duality_gap = semilogy(E, 'ko-','DisplayName','Duality Gap');
hold off

title('Error and Duality Gap')
xlabel('step')
ylabel('error')
legend([absolute_gap, relative_gap, duality_gap], 'Location','best')

saveas(gcf, fullfile('results', date, strcat('error_FW_', line_search, '.png')))

end