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

gap_RE = abs(fx - f_star) / max(1,abs(f_star));
gap_RDG = E ./ max(1,abs(fx(1:end-1)));

gcf = figure('Name', strcat('curve_FW_', line_search));
duality_gap = semilogy(0:length(E)-1, E, 'k-','DisplayName','Duality Gap (DG)');
hold on
relative_duality_gap = semilogy(0:length(E)-1, gap_RDG, 'b-','DisplayName','Relative DG');
relative_error = semilogy(0:length(fx)-1, gap_RE, 'r-','DisplayName','Relative Error');
hold off

title('Error and Duality Gap')
xlabel('step')
ylabel('error')
legend([duality_gap, relative_duality_gap, relative_error], 'Location','northeast')
fontsize(gcf,scale=1.4)

saveas(gcf, fullfile('results', date, strcat('error_FW_', line_search, '.png')))

end