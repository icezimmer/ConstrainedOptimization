function [indices, partition] = PartitionDomain(P)
%{
Give the indices belonging to the simplices with at least two vertices and the partition of these indices.
Input:
    Q  : (matrix) nxn positive semi-definite
    q  : (vector) of length n
    P  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
Output:
    indices   : (vector) indices belongind to the simplices with at least two vertices
    partition : (cell-array) each array represents a simplex
%}

simplices = sum(P,2)>1;

P_no = P(~simplices,:);
indices_no = sum(P_no,1);
indices = ~indices_no;

% Take the non-zero indices (indices in I_k)
P = P(simplices, :);
partition = cell(1,size(P,1));
k=1;
for Pk = P'
    partition{k}=find(Pk)';
    k=k+1;
end

end