function plotLS(Q, q, x, d, alpha, alphaStart)
f = @(t) (x+t*d)'*Q*(x+t*d) + q'*(x+t*d);

t=0;
y(1)=f(t);
delta=alphaStart/100;
i=1;
while(t < alphaStart)
    i=i+1;
    t=t+delta;
    y(i)=f(t);
end

figure('Name','Tomography');
plot(linspace(0,alphaStart,i),y, 'k')
hold on

plot(alpha, f(alpha), 'r*')
hold off