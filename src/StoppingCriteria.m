function condition = StoppingCriteria(fx, f_star, eps_R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%{
if nargin<5
    condition = abs(f(x) - f(x_new)) < eps_R * (1 + abs(f(x)));
elseif nargin==5
    condition = (f(x) - f_star) < eps_R * (1 + abs(f_star));
end
%}

condition = (fx - f_star) < eps_R * (1 + abs(f_star));

end