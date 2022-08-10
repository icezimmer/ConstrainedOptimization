function alpha = StepSizeSelection(Q, q, x, d, i, step_size_method)
%{
Plot the tomography
Input:
    Q                : (matrix) nxn positive semi-definite
    q                : (vector) of length n
    x                : (vector) start point
    d                : (vector) descent direction
    i                : (integer) iteration number
    step_size_method : (string) method for the step size selection
%}

if isequal(step_size_method, 'Exact')
    alpha = ExactLineSearch(Q, q, x, d); 
elseif isequal(step_size_method,'Default')
    alpha = 2/(i + 2);
else
    error("Step Size Method name wrong")
end


end