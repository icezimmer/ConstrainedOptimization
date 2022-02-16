function PrintParameters(n, dim_ker, spectral_radius, K_plus, date)
%{
Print parameters of the instance
Input:
    n               : (integer) dimension of the space
    spectral_radius : (float) spectral radius of the matrix Q
    K_plus          : (integer) number of simplices with at least 2 vertices
    date            : (float) date for saving the results 
%}

path = fullfile('results', date, 'parameters.txt');
fid = fopen(path, 'w+');
fprintf(fid, 'n:%d\ndim_ker:%d\nspectral_radius:%.4f\nK_plus:%d\n', n, dim_ker, spectral_radius, K_plus);
fclose(fid);

end