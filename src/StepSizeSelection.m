function [alpha, alpha_start] = StepSizeSelection(Q, q, x, d, eps_ls, i, step_size_method)
%{
Plot the tomography
Input:
    Q                : (matrix) nxn positive semi-definite
    q                : (vector) of length n
    x                : (vector) start point
    d                : (vector) descent direction
    eps_ls           : (float) stop criterion for the line search
    i                : (integer) iteration number
    step_size_method : (string) method for the step size selection
%}

alpha_start = 1;

if isequal(step_size_method,'LBM')
    alpha_start = StartLineSearch(Q, q, x, d, eps_ls);
    if (alpha_start <= 1)
        alpha = LineSearchLBM(Q, q, x, d, alpha_start, eps_ls);
    else
        alpha = 1;
    end  
elseif isequal(step_size_method,'QBM')
    alpha_start = StartLineSearch(Q, q, x, d, eps_ls);
    if (alpha_start <= 1)
        alpha = LineSearchQBM(Q, q, x, d, alpha_start, eps_ls);
    else
        alpha = 1;
    end 
elseif isequal(step_size_method, 'NM')
    alpha_start = StartLineSearch(Q, q, x, d, eps_ls);
    if (alpha_start <= 1)
        alpha = LineSearchNM(Q, q, x, d, alpha_start, eps_ls);
    else
        alpha = 1;
    end  
elseif isequal(step_size_method,'Default')
    alpha = 2/(i + 2);
end

end