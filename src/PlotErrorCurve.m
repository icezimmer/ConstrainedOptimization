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

gap_A = abs(fx - f_star);
gap_R = gap_A / max(1,abs(f_star));

gcf = figure('Name', strcat('curve_FW_', line_search));
absolute_gap = semilogy(0:length(fx)-1, gap_A, 'b-','DisplayName','Absolute Error');
hold on
relative_gap = semilogy(0:length(fx)-1, gap_R, 'r-','DisplayName','Relative Error');
duality_gap = semilogy(0:length(E)-1, E, 'k-','DisplayName','Duality Gap');
hold off

title('Error and Duality Gap')
xlabel('step')
ylabel('error')
legend([absolute_gap, relative_gap, duality_gap], 'Location','best')

saveas(gcf, fullfile('results', date, strcat('error_FW_', line_search, '.png')))

end