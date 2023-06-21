function SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, num_vertex, date)
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
save(path, "n","dim_ker","spectral_radius","density","K_plus","K_avg","num_vertex")

end