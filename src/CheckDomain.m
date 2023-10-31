function is_in = CheckDomain(x, P)
    %{
        Check if the point x is in the domain, i.e. the point satisfies all K constraints
        Input:
            x     : (vector) point to check
            P     : (matrix) Kxn, each row represents a constraint
        Output:
            is_in : (logical) true iff x satisfies all the K constraints
    %}

    x = x(:);
    [K, ~] = size(P); 
    is_in = (isequal(single(P * x), ones(K, 1, 'like', x)) && all(x >= 0));
end

