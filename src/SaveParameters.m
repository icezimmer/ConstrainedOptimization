function SaveParameters(n, dim_ker, spectral_radius, density, K_plus, date)
%{
Save parameters of the instance
Input:
    n               : (integer) dimension of the space
    spectral_radius : (float) spectral radius of the matrix Q
    density         : (float) density of the matrix Q
    K_plus          : (integer) number of simplices with at least 2 vertices
    date            : (float) date for saving the results 
%}

path = fullfile('results', date, 'parameters.mat');
save(path, "n","dim_ker","spectral_radius","density","K_plus")

end