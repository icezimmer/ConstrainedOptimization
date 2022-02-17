function PlotConvergenceRate(fx, E, line_search, date)
%{
Plot the log-log optimization curve to estimate the convergence rate
Input:
    fx          : (vector) function values sequence
    E           : (vector) error values sequence
    line_search : (string) line search method
    date        : (string) date for saving figures
%}

%addpath symlog_dir

gcf = figure('Name', strcat('rate_FW_', line_search));

tiledlayout(2,1)

nexttile
plot(fx,'go-')
symlog('y')
title('Rate of convergence (primal)')
xlabel('step')
ylabel('f(x)')
% Do nothing else to get just exponents.  Otherwise:
%yt = get(gca,'YTick')';
% Or, for scientific notation
%set(gca,'YTickLabel',strcat('10^',cellstr(num2str(abs(yt)))))

nexttile
semilogy(E, 'mo-')
title('Rate of convergence (dual)')
xlabel('step')
ylabel('duality gap')

saveas(gcf, fullfile('results', date, strcat('rate_FW_', line_search, '.png')))

end

