function condition = StoppingCriteria(fx, f_star, eps_R)
%{
Stopping Criteria for the numerical optimization
Input:
    fx     : (vector) function values sequence
    f_star : (float) optimal value of the function in the domainr
    eps_R  : (float) maximum relative error for stop condition
Output:
    condition : (logical) true if the stop condition is verified
%}

%{
if nargin<5
    condition = abs(f(x) - f(x_new)) < eps_R * (1 + abs(f(x)));
elseif nargin==5
    condition = (f(x) - f_star) < eps_R * (1 + abs(f_star));
end
%}

condition = (fx - f_star) < eps_R * (1 + abs(f_star));

end