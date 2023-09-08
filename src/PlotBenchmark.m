function PlotBenchmark(Histories, f_star, iterative_algorithms, colors, date)
%{
Plot the error curve
Input:
    fx          : (vector) function values sequence
    f_star      : (float) optimal value of the function in the domain
    E           : (vector) error values sequence
    line_search : (string) line search method
    date        : (string) date for saving figures
%}

gcf1 = figure('Name', 'Comparison');
hold on
for i = 1:length(iterative_algorithms)
    history = Histories{i};
    history_f = history.f;
    gap_RE = abs(history_f - f_star) / max(1,abs(f_star));
    semilogy(0:length(gap_RE)-1, gap_RE, 'Color', colors(i),'DisplayName', iterative_algorithms(i));
end
hold off
set(gca, 'YScale', 'log');
title('Primal Comparison')
xlabel('step')
ylabel('relative error')
legend('Location','northeast')
grid on
fontsize(gcf1,scale=1.4)
saveas(gcf1, fullfile('results', date, 'primal_comparison.png'))

gcf2 = figure('Name', 'Comparison');
hold on
for i = 1:length(iterative_algorithms)
    disp(iterative_algorithms(i))
    history = Histories{i};
    history_f = history.f;
    history_dg = history.dg;
    gap_RDG = history_dg ./ max(1,abs(history_f(1:end-1)));
    semilogy(0:length(gap_RDG)-1, gap_RDG, 'Color', colors(i),'DisplayName', iterative_algorithms(i));
end
hold off
set(gca, 'YScale', 'log');
title('Dual Comparison')
xlabel('step')
ylabel('relative error')
legend('Location','northeast')
grid on
fontsize(gcf2,scale=1.4)
saveas(gcf2, fullfile('results', date, 'dual_comparison.png'))
