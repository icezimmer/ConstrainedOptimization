function [x_min_original_dim, f_min, elapsed_time, num_steps, method, variant, primal_error, dual_error, converging, feasible, history] = FrankWolfe(Q, q, P, varargin)
    %{
        FrankWolfe computes the minimum of a quadratic function in a convex compact domain. 
        Input:
            Q                  : (matrix) nxn positive semi-definite
            q                  : (vector) of length n
            P                  : (matrix) Kxn, K is the number of subset I_k and P(k,j) = 1 iff j is in I_k
        Optional input:
            variant            : (string) "Away-step" or "Standard"
            eps_RDG            : (float) maximum dual error for stop condition
            eps_RE             : (float) maximum primal error to consider the solution correct
            max_steps          : (integer) stop criterion (max number of steps for Frank Wolfe)
            tomography         : (logical) plot or not the tomography for each step
            error_plot         : (logical) plot or not the error curve
            date               : (string) date for saving the results
        Output:
            x_min_original_dim : (vector) minimum of the function in the original domain
            f_min              : (vector) min of the function
            elapsed_time       : (float) time elapsed for the computation
            num_steps          : (integer) number of steps for the convergence
            method             : (string) name of the method
            variant            : (string) variant of the method
            primal_error       : (float) relative error
            dual_error         : (float) relative duality gap
            converging         : (logical) method converge or not in the primal perspective
            feasible           : (logical) solution is feasible or not
            history            : (struct) history of the algorithm
    %}

    numvarargs = length(varargin);
    if numvarargs > 8
        error('myfuns:FrankWolfe:TooManyInputs', ...
            'requires at most 8 optional inputs');
    end

    if numvarargs == 8
        f_star = varargin{8};
    elseif numvarargs < 8
        [~, f_star] = Oracle(Q, q(:), P);
    end

    % Set defaults for optional inputs
    optargs = {"Away-step",1e-6,1e-10,10000,false,false,string(datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss'))};
    optargs(1:numvarargs-1) = varargin(1:numvarargs-1);
    [variant, eps_RDG, eps_RE, max_steps, tomography, error_plot, date] = optargs{:};

    if ~exist(fullfile('results',date), 'dir')
        mkdir(fullfile('results',date));
    end

    method = "FW";

    if isequal(variant,'Away-step')
        disp('Away-step Frank-Wolfe algorithm')
    elseif isequal(variant,'Standard')
        disp('Standard Frank-Wolfe algorithm')
    else
        error('Wrong variant name')
    end

    % Dimension of the oringinal space before the semplification
    original_dim = length(q);

    % Semplify the task restricting the optimization on polytopes with at least two vertices
    [Q, q, c, P, partition, fixed] = SemplifyTask(Q, q, P);
    % New function f
    if ~isempty(partition)
        f = @(x) x'*Q*x + q'*x + c;
    else % if K = n
        f = @(x) c;
    end


    % Construct the starting point for the Frank-Wolfe algorithm
    x_start = StartingPoint(P);

    % Force q and x to be column vectors
    q = q(:);
    x = x_start(:);

    tic
    i = 0;
    history_f = f(x);
    duality_gap = NaN;
    history_dg = zeros(0,1);
    % Iterate until convergence
    while (~StoppingCriteria(history_f(end), duality_gap, eps_RDG) && i < max_steps)

        [d, ~, duality_gap] = LinearizationMinimizer(Q, q, x, partition);
        alpha_max = 1;

        if isequal(variant, 'Away-step')
            [d_a, ~, duality_gap_a, alpha_max_a] = LinearizationMaximizer(Q, q, x, partition);
            
            if duality_gap_a > duality_gap
                % Away-step
                d = d_a;
                duality_gap = duality_gap_a;
                alpha_max = alpha_max_a;
            end
        end

        alpha = ExactLineSearch(Q, d, duality_gap, alpha_max);
        
        % Plot tomography
        if tomography
            Tomography(Q, q, x, d, alpha, alpha_max, i, duality_gap)
        end
        
        % Append new value of the duality_gap
        history_dg(i+1) = duality_gap;

        % Upgrade the point
        x = NewPoint(x, d, alpha);

        % Append new value of the funtion
        history_f = cat(2, history_f, f(x));
        
        i = i + 1;
    end

    % Elapsed time for computation
    elapsed_time = toc;

    % Minimum belonging to the semplified domain
    x_min = x;
    % Minimum belonging to the original domain
    x_min_original_dim = ones(original_dim,1);
    x_min_original_dim(~fixed) = x_min;
    % Minimum value of the function
    f_min = history_f(end);

    % Number of steps
    num_steps = i;

    % Feasibility of the solution
    feasible = CheckDomain(x_min, P);

    % Primal precision of the approximation
    [precise, primal_error] = PrimalPrecision(f_min, f_star, eps_RE);

    % Primal convergence
    converging = precise & feasible;

    % Compute the relative duality gap (dual error)
    dual_error = duality_gap / max(1,abs(f_min));

    history.f = history_f;
    history.dg = history_dg;

    if tomography
        disp(['it. ', num2str(i), ', f(x) = ', num2str(f_min)])
    end

    % Plot the error
    if error_plot
        PlotErrorCurve(history, f_star, variant, date)
    end  
end