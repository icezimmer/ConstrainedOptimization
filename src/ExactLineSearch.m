function alpha = ExactLineSearch(Q, d, duality_gap)
%{
Compute the exact line search
Input:
    Q           : (matrix) nxn positive semi-definite
    d           : (vector) descent direction
    duality_gap : (float) opposite of scalar product between the gradient in x and the descent direction
Output:
    alpha : (float) step-size value
%}

% If the denominator is zero alpha will be infinity so we'll fix alpha to 1
alpha = duality_gap / (2*d'* Q *d);

if alpha > 1
    alpha = 1;
end


end