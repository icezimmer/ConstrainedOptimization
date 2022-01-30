function alpha = LineSearch(Q, q, x, d, eps_ls, i, line_search)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if isequal(line_search,'LBM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchLBM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end  
elseif isequal(line_search,'QBM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchQBM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end 
elseif isequal(line_search, 'NM')
    alphaStart = StartLineSearch(Q, q, x, d, eps_ls);
    if (alphaStart <= 1)
        alpha = LineSearchNM(Q, q, x, d, alphaStart, eps_ls);
    else
        alpha = 1;
    end  
elseif isequal(line_search,'Default')
    alpha = 2/(i + 2);
end

end

