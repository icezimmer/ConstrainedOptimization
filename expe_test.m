%{
    Perform Frank-Wolfe to several instances of the task
%}

addpath src

force_non_point_simplices = true;
variant = "Away-step";
tomography = false;
error_plot = true;
seed = 1;
eps_RDG = 1e-6; max_steps = 1e6;
eps_RE = 1e-10;

% Lists of values for the parameters
n_list = {100};
K_list = {@(n)1,@(n)floor(0.1*n),@(n)floor(n*exp(-1)),@(n)floor(0.5*n)};
actv_list = {0,0.5,1};
dim_ker_list = {@(n)0.1*n};
spectral_radius_list = {10};
lambda_min_list = {1};
density_list = {1};

num_trials=length(n_list)*length(K_list)*length(actv_list)*length(dim_ker_list)*length(spectral_radius_list)*length(lambda_min_list)*length(density_list);
trial=0;
for i_=1:length(n_list)
    n = n_list{i_};
    for j_=1:length(K_list)
        K = K_list{j_};
        K = K(n);
        for k_=1:length(actv_list)
            actv = actv_list{k_};
            for l_=1:length(dim_ker_list)
                dim_ker=dim_ker_list{l_};
                dim_ker=dim_ker(n);
                for m_=1:length(spectral_radius_list)
                    spectral_radius = spectral_radius_list{m_};
                    for n_=1:length(lambda_min_list)
                        lambda_min=lambda_min_list{n_};
                        for o_=1:length(density_list)
                            density=density_list{o_};
                            pause(2)
                            [Q, q, P, K_plus, K_avg, num_vertex, norm_q, date] = GenerateInstance(n, K, force_non_point_simplices, actv, dim_ker, spectral_radius, lambda_min, density, seed);
                            SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_ker, spectral_radius, lambda_min, density, norm_q, seed, date)
                            SaveMatrices(Q, q, P, date)
                            [x_min, f_min, elapsed_time, num_steps, method, variant, err, converging, feasible, duality_gap, history] = FrankWolfe(Q, q, P, variant, eps_RDG, eps_RE, max_steps, tomography, error_plot, date);
                            SaveTestResults(x_min, f_min, elapsed_time, num_steps, method, variant, err, converging, feasible, duality_gap, history, date)
                            trial = trial+1;
                            disp([num2str(trial*100/num_trials),'%'])
                        end
                    end
                end
            end
        end
    end
end
