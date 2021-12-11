function [x_min, f_min] = GlobalMinimum(Q, q)
%{
Input:
    Q : (matrix) nxn positive semi-definite
    q : (vector) of length n
Output:
    x_min : (vector) global argmin of the function
    f_min : (vector) global minimum of the function
%}

% Function f
f = @(x) x'*Q*x + q'*x;

if (rank(Q) == rank([Q, q]))
    x_min = Q \ (-q);
    f_min = f(x_min);
end

end

