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

% Subvector of z that violates the constraints
% Subvector belongs to [-1, 0]^|I_k| or [1, 2]^|I_k|
negative = rand(K_out,1) <= 0.5;
center_negative = -1/2;
center_positive =  3/2;
radius = 1/2;
center = center_negative*negative + center_positive*(~negative);
for k = 1 : K_out
    indices_k = find(P(k,:));
    x_k = center(k)-radius + 2*radius*rand(length(indices_k),1);
    z(indices_k) = x_k;
end 

% Subvector of z that satisfies the constraints
for k = K_out+1 : K
    indices = find(P(k,:));
    x_k = randfixedsum(length(indices),1,1,0,1);
    z(indices) = x_k;
end

q = 2*(-z'*Q)';

end
