function PlotBenchmark(Histories, f_star, iterative_algorithms, dual_comparison, colors, date)
%{
Plot the error curve
Input:
    fx          : (vector) function values sequence
    f_star      : (float) optimal value of the function in the domain
    E           : (vector) error values sequence
    line_search : (string) line search method
    date        : (string) date for saving figures
%}

gcf1 = figure('Name', 'Primal Comparison');

hold on
for i = 1:length(iterative_algorithms)
    history = Histories{i};
    history_f = history.f;
    gap_RE = cummin(abs(history_f - f_star) / max(1,abs(f_star)));
    semilogy(0:length(gap_RE)-1, gap_RE, 'Color', colors(i),'DisplayName', iterative_algorithms(i), 'LineWidth',2);
end
hold off
set(gca, 'YScale', 'log');
title('Primal Comparison')
xlabel('step')
ylabel('primal error')
legend('Location','northeast')
fontsize(gcf1,scale=1.4)
grid on
% Get the figure's position to calculate the width and height
figPosition = get(gcf1, 'Position');
width = figPosition(3);
height = figPosition(4);
% Set the paper size to match the figure size
set(gcf1, 'PaperUnits', 'points');
set(gcf1, 'PaperSize', [width, height]);
set(gcf1, 'PaperPosition', [0, 0, width, height]);
exportgraphics(gcf1, fullfile('results', date, 'primal_comparison.pdf'), 'Resolution', 300);

if dual_comparison
    gcf2 = figure('Name', 'Dual Comparison');
    hold on
    for i = 1:length(iterative_algorithms)
        history = Histories{i};
        history_f = history.f;
        history_dg = history.dg;
        gap_RDG = cummin(history_dg ./ max(1,abs(history_f(1:end-1))));
        semilogy(0:length(gap_RDG)-1, gap_RDG, 'Color', colors(i),'DisplayName', iterative_algorithms(i), 'LineWidth',2);
    end
    hold off
    set(gca, 'YScale', 'log');
    title('Dual Comparison')
    xlabel('step')
    ylabel('dual error')
    legend('Location','northeast')
    fontsize(gcf2,scale=1.4)
    grid on
    % Get the figure's position to calculate the width and height
    figPosition = get(gcf2, 'Position');
    width = figPosition(3);
    height = figPosition(4);
    % Set the paper size to match the figure size
    set(gcf2, 'PaperUnits', 'points');
    set(gcf2, 'PaperSize', [width, height]);
    set(gcf2, 'PaperPosition', [0, 0, width, height]);
    exportgraphics(gcf2, fullfile('results', date, 'dual_comparison.pdf'), 'Resolution', 300);
end
