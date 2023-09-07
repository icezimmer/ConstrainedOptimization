function [converging, err] = ConvergingError(f_x, f_star, eps_E, varargin)
%{
Stopping Criteria for the numerical optimization
Input:
    f_x     : (float) function value at point x
    f_star  : (float) optimal value of the function in the domain
    eps_E  : (float) maximum error for converging
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
    err = (f_x - f_star) / max(1,abs(f_star));
elseif isequal(type,'Absolute')
    err = (f_x - f_star) / max(1,abs(f_star));
else
    error("Wrong stopping criteria")
end

converging = err < eps_E;

end