function x_new = NewPoint(x_old, d, alpha)
%{
Compute the new point given the old point, the descent direction and the step size.
Input:
    x_old : (vector) old point
    d     : (vector) descent direction
    alpha : (float) step size
Output:
    x_new : (vector) new point
%}
x_new = x_old + alpha * d;
end

