function SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_Ker, spectral_radius, lambda_min, density, norm_q, seed, date)
    %{
        Save the parameters of the instance in a .mat file
        Input:
            n               : (integer) dimension of the space
            K_plus          : (integer) number of simplices with at least 2 vertices
            K_avg           : (integer) average of dimensions of the simplices with at least 2 vertices
            num_vertex      : (integer) number of vertices
            actv            : (integer) fraction of active constraints
            dim_ker         : (integer) ker dimension of the matrix Q
            spectral_radius : (float) spectral radius of the matrix Q
            lambda_min      : (float) min strictly positive eigenvalue of the matrix Q
            density         : (float) density of the matrix Q
            norm_q          : (float) norm of the vector q
            seed            : (integer) seed for the random generator
            date            : (float) date for saving the results 
    %}

    path = fullfile('results', date, 'parameters.mat');
    if dim_Ker == n-1
        lambda_min = spectral_radius;
    end
    parameters.n=n;
    parameters.K_plus=K_plus;
    parameters.K_avg=K_avg; 
    parameters.num_vertex=num_vertex; 
    parameters.actv=actv; 
    parameters.dim_Ker=dim_Ker;
    parameters.spectral_radius=spectral_radius; 
    parameters.lambda_min=lambda_min; 
    parameters.density=density;
    parameters.norm_q = norm_q;
    parameters.seed = seed;

    save(path, "parameters")
end