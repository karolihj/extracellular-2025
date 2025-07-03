function P = set_up_param_from_file(param, G)
% P = set_up_param_from_file(param, G)

P = param*ones(1, G.Nm);

if isfield(G, 'parameters_from_file')
    for n=1:length(G.parameters_from_file)
        param_idx = find(ismember(G.param_names, G.parameters_from_file{n}));
        param_values = importdata(G.parameter_files{n});
        P(param_idx,:) = param_values(1:G.Nm);
    end
end

end

