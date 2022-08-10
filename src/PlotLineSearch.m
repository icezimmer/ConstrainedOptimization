function PlotLineSearch(Q, q, x, d, alpha, i)
%{
Plot the tomography
Input:
    Q     : (matrix) nxn positive semi-definite
    q     : (vector) of length n
    x     : (vector) starting point
    d     : (vector) descent direction
    alpha : (float) coefficient for the descent direction
    i     : (integer) iteration number
%}

f = @(t) (x+t*d)'*Q*(x+t*d) + q'*(x+t*d);

t=0;
y=f(t);
delta=(2*alpha)/100;
stop = (2*alpha)*(2*alpha<1) + 1*(2*alpha >=1);
while(t < stop)
    t=t+delta;
    y=cat(1,y,f(t));
end

figure('Name','Tomography');
plot(linspace(0,stop,length(y)),y, 'k')
title(['From ', 'f(x(', num2str(i), '))', '  to  ', 'f(x(', num2str(i+1), '))'])
xlabel('alpha')
ylabel('f(x)')
hold on

plot(alpha, f(alpha), 'r*')
hold off

end