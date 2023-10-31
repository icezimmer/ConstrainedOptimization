function alpha = ExactLineSearch(Q, d, duality_gap, alpha_max)
    %{
        Compute the exact line search
        Input:
            Q           : (matrix) nxn positive semi-definite
            d           : (vector) descent direction
            duality_gap : (float) opposite of scalar product between the gradient in x and the descent direction
            alpha_max   : (float) maximum step-size value
        Output:
            alpha       : (float) step-size value
    %}

    alpha = min(alpha_max, duality_gap / (2*d'* Q *d));
end