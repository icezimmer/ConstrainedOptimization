function plotMOMENTUM(Q, q, x, d_old, par_momentum, d_new, alpha_new)
%{
Plot the tomography (momentum version)
Input:
    Q            : (matrix) nxn positive semi-definite
    q            : (vector) of length n
    x            : (vector) start point
    d_old        : (vector) old direction
    par_momentum : (float) momentum coefficient
    d_new        : (vector) new direction
    alpha_new    : (float) coefficient for the new direction
%}

d = (alpha_new * d_new) + (par_momentum * d_old); 

f = @(t) (x+t*d)'*Q*(x+t*d) + q'*(x+t*d);

t=0;
y(1)=f(t);
delta=1/100;
i=1;
while(t < 1)
    i=i+1;
    t=t+delta;
    y(i)=f(t);
end

figure('Name','Tomography with momentum');
plot(linspace(0,1,i),y, 'k')
hold on

plot(1, f(1), 'r*')
hold off