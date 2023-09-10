function PlotErrorCurve(history, f_star, variant, date)
%{
Plot the error curve
Input:
    fx      : (vector) function values sequence
    f_star  : (float) optimal value of the function in the domain
    DG      : (vector) duality gap values sequence
    variant : (string) variant of the FW algorithm
    date    : (string) date for saving figures
%}

f_x = history.f;
DG = history.dg;

gap_RE = abs(f_x - f_star) / max(1,abs(f_star));
gap_RDG = cummin(DG ./ max(1,abs(f_x(1:end-1))));

gcf = figure('Name', strcat('curve_FW_', variant));
relative_duality_gap = semilogy(0:length(DG)-1, gap_RDG, 'Color',"#0072BD", 'DisplayName','Relative Duality Gap', 'LineWidth',2);
hold on
relative_error = semilogy(0:length(f_x)-1, gap_RE, 'Color', "#D95319",'DisplayName','Relative Error', 'LineWidth',2);
hold off

title('Primal and Dual Error')
xlabel('step')
ylabel('error')
legend([relative_duality_gap, relative_error], 'Location','northeast')
fontsize(gcf,scale=1.4)
grid on

% Get the figure's position to calculate the width and height
figPosition = get(gcf, 'Position');
width = figPosition(3);
height = figPosition(4);
% Set the paper size to match the figure size
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [width, height]);
set(gcf, 'PaperPosition', [0, 0, width, height]);

% Save the figure
print(gcf, fullfile('results', date, strcat('error_FW_', variant, '.pdf')), '-dpdf');

end