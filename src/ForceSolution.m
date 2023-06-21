function [z, q] = ForceSolution(Q, P, actv)
% - actv (real, scalar, default 0.5): how many box constraints (as a
%   fraction of the number of variables n of the problems) the
%   unconstrained optimum will violate, and therefore we expect to be
%   active in the constrained optimum; note that there is no guarantee that
%   exactly acvt constraints will be active, they may be less or (more
%   likely) more, except when actv = 0 because then the unconstrained
%   optimum is surely feasible and therefore it will be the constrained
%   optimum as well

if actv<0 || actv>1
    error("The fracion of active constraints must be in [0,1]")
end

[K,n] = size(P);
z = zeros(n,1);

% The unconstrained optimum z violates the first floor(actv*K) simplices (the respective x_k are in boundary of P_k)
K_out = floor(actv * K);

% subvector of z that violates the constraints
for k = 1 : K_out
    indices = find(P(k,:));
    % x_k = zeros(length(indices),1);
    % negative = rand(length(indices),1) <= 0.5;
    % x_k(negative) = -10 * rand(nnz(negative),1);
    % x_k(~negative) = 1 + 10 * rand(nnz(~negative),1);
    center = -2;
    radius = 1;
    x_k = center-radius + 2*radius*rand(length(indices),1);
    z(indices) = x_k;
end 

% subvector of z that satisfies the constraints
for k = K_out+1 : K
    indices = find(P(k,:));
    x_k = randfixedsum(length(indices),1,1,0,1);
    z(indices) = x_k;
end

q = 2*(-z'*Q)';

end
