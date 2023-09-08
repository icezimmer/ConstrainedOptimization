function PlotErrorCurve(f_x, f_star, DG, variant, date)
%{
Plot the error curve
Input:
    fx      : (vector) function values sequence
    f_star  : (float) optimal value of the function in the domain
    DG      : (vector) duality gap values sequence
    variant : (string) variant of the FW algorithm
    date    : (string) date for saving figures
%}

gap_RE = abs(f_x - f_star) / max(1,abs(f_star));
gap_RDG = DG ./ max(1,abs(f_x(1:end-1)));

gcf = figure('Name', strcat('curve_FW_', variant));
relative_duality_gap = semilogy(0:length(DG)-1, gap_RDG, 'b-','DisplayName','Relative Duality Gap');
hold on
relative_error = semilogy(0:length(f_x)-1, gap_RE, 'r-','DisplayName','Relative Error');
hold off

title('Primal and Dual Error')
xlabel('step')
ylabel('error')
legend([relative_duality_gap, relative_error], 'Location','northeast')
fontsize(gcf,scale=1.4)

saveas(gcf, fullfile('results', date, strcat('error_FW_', variant, '.png')))

end