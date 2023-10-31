function [Q_new, q_new, c, P_new, partition, fixed] = SemplifyTask(Q, q, P)
    %{
        Semplify the task by removing the fixed variables.
        From minimizing f(x)=x'*Q*x+q'*x in R^n to minimizing f_new(y)=y'*Q_new*y+q_new'*y+c in R^m,
        where m is the number of free variables.
        Input:
            Q         : (matrix) nxn positive semi-definite
            q         : (vector) of length n
            P         : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
        Output:
            Q_new     : (matrix) nxn positive semi-definite
            q_new     : (vector) of length n
            c         : (float) constant
            P_new     : (matrix) Kxn,
            partition : (cell) of length K, partition{k} is the subset I_k
            fixed     : (vector) of length n, fixed(j) = 1 iff j is fixed
    %}
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
end