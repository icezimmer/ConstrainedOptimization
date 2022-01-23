function plotCURVE(fx, E, line_search, beta)
%{
Plot the tomography
Input:
    fx          : (vector) function values sequence
    E           : (vector) error values sequence
    line_search : (string) line search method
    beta        : (float) momentum coefficient
%}

figure('Name', strcat(line_search, " with momentum = ", num2str(beta)));

tiledlayout(2,1)

nexttile
plot(fx, 'bo-')
title('Optimization curve')
xlabel('step')
ylabel('f(x)')

nexttile
plot(abs(E), 'ro-')
title('Error curve')
xlabel('step')
ylabel('error')

end