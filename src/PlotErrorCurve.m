function PlotErrorCurve(history, f_star, variant, date)
    %{
        Plot the error curve
        Input:
            history : (struct) history of the algorithm
            f_star  : (float) optimal value of the function in the domain
            variant : (string) variant of the FW algorithm
            date    : (string) date for saving figures
    %}

    % Get the function values and the duality gap
    f_x = history.f;
    DG = history.dg;

    % Compute the primal and dual error
    gap_RE = cummin(abs(f_x - f_star) / max(1,abs(f_star)));
    gap_RDG = cummin(DG ./ max(1,abs(f_x(1:end-1))));

    gcf = figure('Name', strcat('curve_FW_', variant));
    relative_duality_gap = semilogy(0:length(DG)-1, gap_RDG, 'Color',"#0072BD", 'DisplayName','Dual Error', 'LineWidth',2);
    hold on
    relative_error = semilogy(0:length(f_x)-1, gap_RE, 'Color', "#D95319",'DisplayName','Primal Error', 'LineWidth',2);
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
    exportgraphics(gcf, fullfile('results', date, strcat('error_FW_', variant, '.pdf')), 'Resolution', 300);
end