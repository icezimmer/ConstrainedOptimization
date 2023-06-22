function [Q_new, q_new, c, P_new, partition, fixed] = SemplifyTask(Q, q, P)
simplices = sum(P,2)>1;

free = logical(sum(P(simplices,:),1))';
fixed = ~free;

% Save in partition the subsets I_k with at least two indices
P_new = P(simplices, free);
partition = cell(1,size(P_new,1));
k=1;
for Pk = P_new'
    partition{k}=find(Pk)';
    k=k+1;
end

Q_new = Q(free,free);
q_new = q(free) + 2*(sum(Q(free,fixed),2));
c = sum(Q(fixed,fixed), 'all')+sum(q(fixed));