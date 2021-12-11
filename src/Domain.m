function is_in = Domain(x, P)
%{
Check if the point x is in the domain, i.e. the point satisfies all K constraints
Input:
    x : (vector) point to check
    P : (matrix) Kxn, each row represent a constraint
Output:
    is_in : (logical) true iff x satisfies all the K constraints
%}

x = x(:);
[K, ~] = size(P); 
is_in = all(single(P * x) == single(ones(length(K), 1, 'like', x)));

end

