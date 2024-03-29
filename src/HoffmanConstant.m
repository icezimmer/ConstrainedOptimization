function HC = HoffmanConstant(A)
    %{
        Compute the Hoffman Constant of a matrix A
        Input:
            A: (matrix) the matrix
        Output:
            HC: (float) the Hoffman Constant
    %}
    disp("Compute the Hoffman Constant")
    HC = -inf;
    for k=1:rank(full(A))
        BB = FindLinearlyIndependentRows(full(A),k);
        for j=1:length(BB)
            B = BB{j};
            new_value = 1/eigs(B*B',1,'smallestabs');
            if new_value > HC
                HC = new_value;
            end
        end
    end
end
