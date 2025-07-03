function [V_m, U_m] = solve_system(G, mesh)
% [V_m, U_m] = solve_system(G, mesh)

% Set up matrix
A = set_up_matrix(G, mesh);

% Set up default initial conditions and parameters
states = G.init_states(G.param_filename);
param = G.init_parameters(G.param_filename);

% Set up intial conditions for the entire domain
states = states*ones(1, G.Nm);
t = 0;
V = states(G.V_idx,:)';
U = zeros(G.N, 1);

% Set up cell model parameters
P = set_up_stim_param(param, G, mesh);
P = introduce_variation(P, G);

% Set up zero vector
d = zeros(G.N, 1);

% Prepare for saving solutions
if isfield(G, 'DT')
    G.num_save = round(G.Tstop/G.DT) + 1;
    G.save_step = round(G.DT/G.dt);
else
    G.num_save = G.Nt+1;
    G.save_step = 1;
    G.DT = G.dt;
end


% Matrix for saving the solution
V_m = zeros(G.Nm, G.num_save);
V_m(:,1) = V;
U_m = zeros(G.N, G.num_save);
U_m(:,1) = U;

% Run simulation
fprintf('Starting simulation...\n')
print_step = 50*G.dt;
n_print = max(round(print_step/G.dt), 1);
t1 = tic;
t_mem = 0;
t_i = 0;
for n = 1:G.Nt
    
    % Step 1: Solve ode-system for the single cell model (forward Euler)
    for k=1:G.nt
        states = states + G.dt_ode*G.rhs(t, states, P);
        t = t + G.dt_ode;
    end
    V = states(G.V_idx, :)';
    if any(isnan(V)) || any(isinf(V))
        fprintf('ODE solution is inf or nan.\n')
        break;
    end
    
    
    % Step 2: Solve KNM-system 
    b = [V; d];
    X = A\b;
    
    % Extract solution
    V = X(1:G.Nm);
    U = X(G.Nm+(1:G.N));
    
    % Save solution
    if rem(n, G.save_step) == 0
        V_m(:, round(n/G.save_step)+1) = V;
        U_m(:, round(n/G.save_step)+1) = U;
    end
   
    % Update the membrane potential in the single cell model
    states(G.V_idx, :) = V';
    
    % Estimate remaining simulation time
    if rem(n, n_print) == 0
        % Print current point in time
        fprintf('t = %.2f ms. ', t);

        % Print estimated simulation time
        t2 = toc(t1);                % Time usage for n_print time steps
        t_rem = t2*(G.Nt-n)/n_print; % Estimated remaining simulation time
        fprintf('Estimated remaining simulation time: ');
        if t_rem > 86400*2
            fprintf('%.1f days \n', t_rem/86400);
        elseif t_rem > 3600
            fprintf('%.1f h \n', t_rem/3600);
        elseif t_rem > 60
            fprintf('%.1f min \n', t_rem/60);
        else
            fprintf('%.1f sec \n', t_rem);
        end
        t1 = tic;
        
    end
end
fprintf('Simulation done\n\n')

V_new = nan*ones(G.N, G.num_save);
V_new(G.with_myocyte, :) = V_m;
V_m = V_new;


end


