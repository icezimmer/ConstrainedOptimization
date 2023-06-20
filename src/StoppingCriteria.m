function condition = StoppingCriteria(fx, f_star, eps_R, varargin)
%{
Stopping Criteria for the numerical optimization
Input:
    fx     : (float) function value at point x
    f_star : (float) optimal value of the function in the domain
    eps_R  : (float) maximum relative error for stop condition
Output:
    condition : (logical) true if the stop condition is verified
%}

numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:StoppingCriteria:TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {'Relative'};
optargs(1:numvarargs) = varargin;
[type] = optargs{:};

if isequal(type,'Relative')
    condition = (fx - f_star) < eps_R * max(1,abs(f_star));
elseif isequal(type,'Absolute')
    condition = (fx - f_star) < eps_R;
else
    error("Wrong stopping criteria")
end

end