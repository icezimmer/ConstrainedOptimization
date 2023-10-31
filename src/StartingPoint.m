function x_start = StartingPoint(P)
    %{
        Compute the starting point for the Frank-Wolfe algorithm
        Input:
            P : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    %}

    [K, n] = size(P);
    x_start = zeros(n, 1);
    for k = 1 : K
        i = find(P(k,:), 1);
        x_start(i) = 1;
    end
end