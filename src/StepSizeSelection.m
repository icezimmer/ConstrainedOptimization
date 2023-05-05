function alpha = StepSizeSelection(Q, d, duality_gap, alpha_max, i, step_size_method)
%{
Step size selection
Input:
    Q                : (matrix) nxn positive semi-definite
    d                : (vector) descent direction
    duality_gap      : (float) opposite of scalar product between the gradient in x and the descent direction
    i                : (integer) iteration number
    step_size_method : (string) method for the step size selection
%}

if isequal(step_size_method, 'Exact')
    alpha = ExactLineSearch(Q, d, duality_gap, alpha_max); 
elseif isequal(step_size_method,'Standard')
    alpha = 2/(i + 2);
else
    error("Step Size Method name wrong")
end


end