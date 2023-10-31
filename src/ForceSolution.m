function [z, q] = ForceSolution(Q, P, actv)
    %{
        Generates a feasible solution z for the problem
        Input:
            Q    : (matrix) nxn positive semi-definite
            P    : (matrix) Kxn, each row represents a constraint, s.t.
                the domain is D = {x in R^n : P*x = 1, x >= 0}
            actv : (real) fraction of active constraint, i.e. how many constraint the unconstrained optimum will violate, 
                and therefore we expect to be active in the constrained optimum
        Output:
            z : (vector) nx1, unconstrained optimum
            q : (vector) nx1, vector q s.t. f(x) = x'*Q*x + q'*x has the unconstrained optimum z
    %}

    if actv<0 || actv>1
        error("The fracion of active constraints must be in [0,1]")
    end

    [K,n] = size(P);
    z = zeros(n,1);

    % The unconstrained optimum z violates the first floor(actv*K) simplices (the respective x_k are in the 
    % relative boundary of P_k)
    K_out = floor(actv * K);

    % Subvector of z that violates the constraints,
    % in particular the subvector belongs to [-1, 0]^|I_k| or [1, 2]^|I_k|
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

    % Each k-th subvector of z forced to belong to relint(P_k)
    for k = K_out+1 : K
        indices = find(P(k,:));
        relint = false;
        while(~relint)
            x_k = randfixedsum(length(indices),1,1,0,1);
            if nnz(x_k)==length(indices)
                relint = true;
            end
        end
        z(indices) = x_k;
    end

    q = 2*(-z'*Q)';
end
