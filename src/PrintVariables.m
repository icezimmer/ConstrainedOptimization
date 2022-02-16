function PrintVariables(Q, q, P, date)
%{
Print the variables Q, q and P that represents the domain
Input:
    Q    : (matrix) nxn positive semi-definite
    q    : (vector) of length n
    P    : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
    date : (float) date for saving the results 
%}

path = fullfile('results', date, 'variables.txt');

fid = fopen(path, 'w+');

fprintf(fid, 'Q,q,P\n');  % header

dlmwrite(path, Q, 'delimiter', '\t', 'precision', 4, '-append', 'roffset', 1)

dlmwrite(path, q, 'delimiter', '\t', 'precision', 4, '-append', 'roffset', 1)

dlmwrite(path, P, 'delimiter', '\t', 'precision', 4, '-append', 'roffset', 1)

end