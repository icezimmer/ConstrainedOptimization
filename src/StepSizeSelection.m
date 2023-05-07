function alpha = StepSizeSelection(Q, d, duality_gap, alpha_max, i, variant)
%{
Step size selection
Input:
    Q                : (matrix) nxn positive semi-definite
    d                : (vector) descent direction
    duality_gap      : (float) opposite of scalar product between the gradient in x and the descent direction
    i                : (integer) iteration number
    step_size_method : (string) method for the step size selection
%}

if isequal(variant, 'Exact') || isequal(variant, 'Away-step')
    alpha = ExactLineSearch(Q, d, duality_gap, alpha_max); 
elseif isequal(variant,'Standard')
    alpha = 2/(i + 2);
else
    error("Variant name wrong")
end


end