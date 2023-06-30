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

gcf = figure('Name', 'Comparison');
hold on
for i = 1:length(iterative_algorithms)
    gap = abs(Histories{i} - f_star) / max(1,abs(f_star));
%     if isequal(iterative_algorithms(i), 'FW-Away-step')
%         r = 1 / sqrt(2); % max rate of linear convergence for AFW
%         semilogy(0:length(gap)-1, gap(1) * r.^(0:length(gap)-1), 'Color', 'm', 'DisplayName', "Max-Rate-AFW");
%     end
    semilogy(0:length(gap)-1, gap, 'Color', colors(i),'DisplayName', iterative_algorithms(i));
end
hold off
set(gca, 'YScale', 'log');
title('Comparison')
xlabel('step')
ylabel('error')
legend('Location','northeast')
fontsize(gcf,scale=1.4)
saveas(gcf, fullfile('results', date, 'comparison.png'))
