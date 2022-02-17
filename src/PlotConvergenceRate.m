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
% Symlog beacuse fx is not a positive vector
symlog('y')
title('Rate of convergence (primal)')
xlabel('step')
ylabel('f(x)')

nexttile
semilogy(E, 'mo-')
title('Rate of convergence (dual)')
xlabel('step')
ylabel('duality gap')

saveas(gcf, fullfile('results', date, strcat('rate_FW_', line_search, '.png')))

end

