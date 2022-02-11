function PlotLineSearch(Q, q, x, d, alpha, alpha_start, k)
%{
Plot the tomography
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    x           : (vector) start point
    d           : (vector) descent direction
    alpha       : (float) coefficient for the descent direction
    alpha_start : (float) start value for alpha
    k           : (integer) number of iteration
%}

f = @(t) (x+t*d)'*Q*(x+t*d) + q'*(x+t*d);

t=0;
y(1)=f(t);
delta=alpha_start/100;
i=1;
while(t < alpha_start)
    i=i+1;
    t=t+delta;
    y(i)=f(t);
end

figure('Name','Tomography');
plot(linspace(0,alpha_start,i),y, 'k')
title(['From ', 'f(x(', num2str(k), '))', '  to  ', 'f(x(', num2str(k+1), '))'])
xlabel('alpha')
ylabel('f(x)')
hold on

plot(alpha, f(alpha), 'r*')
hold off

end