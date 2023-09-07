function condition = StoppingCriteria(f_x, duality_gap, eps_DG, varargin)
%{
Stopping Criteria for the numerical optimization
Input:
    f_x         : (float) function value at point x
    duality_gap : (float) duality gap at point x
    eps_RDG     : (float) maximum relative duality gap for stop condition
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
    condition = duality_gap < eps_DG * max(1,abs(f_x));
elseif isequal(type,'Absolute')
    condition = duality_gap < eps_DG;
else
    error("Wrong stopping criteria")
end

end