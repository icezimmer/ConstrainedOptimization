function [indices, partition] = PartitionDomain(P)
%{
Give the indices belonging to the simplices with only one vertex and
the partition of the indices belonging to the simplices with at least 2 vertices.
Input:
    Q  : (matrix) nxn positive semi-definite
    q  : (vector) of length n
    P  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
Output:
    indices   : (vector) indices belonging to the simplices with only one vertex
    partition : (cell-array) each array represents a simplex with at least two vertices
%}

simplices = sum(P,2)>1;

indices = logical(sum(P(~simplices,:),1))';

% Save in partition the subsets I_k with at least two indices
P = P(simplices, :);
partition = cell(1,size(P,1));
k=1;
for Pk = P'
    partition{k}=find(Pk)';
    k=k+1;
end

end