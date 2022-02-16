function PlotLineSearch(Q, q, x, d, alpha, alpha_start, i)
%{
Plot the tomography
Input:
    Q           : (matrix) nxn positive semi-definite
    q           : (vector) of length n
    x           : (vector) starting point
    d           : (vector) descent direction
    alpha       : (float) coefficient for the descent direction
    alpha_start : (float) start value for alpha
    i           : (integer) iteration number
%}

f = @(t) (x+t*d)'*Q*(x+t*d) + q'*(x+t*d);

t=0;
y(1)=f(t);
delta=alpha_start/100;
k=1;
while(t < alpha_start)
    k=k+1;
    t=t+delta;
    y(k)=f(t);
end

figure('Name','Tomography');
plot(linspace(0,alpha_start,i),y, 'k')
title(['From ', 'f(x(', num2str(i), '))', '  to  ', 'f(x(', num2str(i+1), '))'])
xlabel('alpha')
ylabel('f(x)')
hold on

plot(alpha, f(alpha), 'r*')
hold off

end