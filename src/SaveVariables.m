function SaveVariables(Q, q, P, date)
%{
Save the variables Q, q and P that represents the domain
Input:
    Q    : (matrix) nxn positive semi-definite
    q    : (vector) of length n
    P    : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    date : (float) date for saving the results 
%}

path = fullfile('results', date, 'variables.mat');
save(path, "Q","q","P");

end