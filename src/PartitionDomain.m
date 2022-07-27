function [indices, partition] = PartitionDomain(P)
%{
Eliminate the siplices of dimension 1 and their respective dimensions.
Input:
    Q  : (matrix) nxn positive semi-definite
    q  : (vector) of length n
    P  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
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