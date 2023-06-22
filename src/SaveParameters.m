function SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_Ker, spectral_radius, lambda_min, density, date)
%{
Save parameters of the instance
Input:
    n               : (integer) dimension of the space
    dim_ker         : (integer) ker dimension of the matrix Q
    spectral_radius : (float) spectral radius of the matrix Q
    density         : (float) density of the matrix Q
    K_plus          : (integer) number of simplices with at least 2 vertices
    K_avg           : (integer) average of dimensions of the simplices with at least 2 vertices
    date            : (float) date for saving the results 
%}

path = fullfile('results', date, 'parameters.mat');
parameters.n=n;
parameters.K_plus=K_plus;
parameters.K_avg=K_avg; 
parameters.num_vertex=num_vertex; 
parameters.actv=actv; 
parameters.dim_Ker=dim_Ker;
parameters.spectral_radius=spectral_radius; 
parameters.lambda_min=lambda_min; 
parameters.density=density;

save(path, "parameters")

end