function [precise, err] = PrimalPrecision(f_x, f_star, eps_E, varargin)
    %{
        Evaluate the primal error (relative error) of the approximation
        Input:
            f_x       : (float) function value at point x
            f_star    : (float) optimal value of the function in the domain
            eps_E     : (float) maximum error for good primal precision
        Optional input:
            type      : (string) type of error (relative or absolute)
        Output:
            condition : (logical) true if the stop condition is verified
    %}

    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfuns:PrimalPrecision:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end

    % Set defaults for optional inputs
    optargs = {'Relative'};
    optargs(1:numvarargs) = varargin;
    [type] = optargs{:};

    if isequal(type,'Relative')
        err = abs(f_x - f_star) / max(1,abs(f_star));
    elseif isequal(type,'Absolute')
        err = abs(f_x - f_star) / max(1,abs(f_star));
    else
        error("Wrong criteria")
    end

    precise = err < eps_E;
end