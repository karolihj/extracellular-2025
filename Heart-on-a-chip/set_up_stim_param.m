function P = set_up_stim_param(param, G, mesh)
% P = set_up_stim_param(param, G, mesh)

[~, stim_idx] = ismember('stim_amplitude', G.param_names);
[~, stim_period_idx] = ismember('stim_period', G.param_names);
param(stim_period_idx) = 1000;
P = param*ones(1, G.Nm);
p_tmp = zeros(G.N, 1);
p_tmp(mesh.to_stim) = G.stim_amp;
p_tmp = p_tmp(G.with_myocyte);
P(stim_idx, :) = p_tmp;

% Stim duration
if isfield(G,  'stim_time')
   [~, stim_time_idx] = ismember('stim_duration', G.param_names); 
   P(stim_time_idx, :) = G.stim_time; 
end

end

