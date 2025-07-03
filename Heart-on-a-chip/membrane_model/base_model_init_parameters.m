function [parameters, varargout] = base_model_init_parameters(filename)
  % % Default parameter values for ODE model: base_model

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  load(filename, 'parameters')

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(86, 1);

    % --- I_Na ---
    parameter_names{1} = 'g_Na';

    % --- I_NaL ---
    parameter_names{2} = 'g_NaL';
    parameter_names{3} = 'thL';

    % --- I_NaK ---
    parameter_names{4} = 'KmKo';
    parameter_names{5} = 'KmNaip';
    parameter_names{6} = 'g_NaK';

    % --- I_Kr ---
    parameter_names{7} = 'g_Kr';

    % --- I_Ks ---
    parameter_names{8} = 'g_Ks';

    % --- I_to ---
    parameter_names{9} = 'g_to';

    % --- I_K1 ---
    parameter_names{10} = 'g_K1';

    % --- I_bCl ---
    parameter_names{11} = 'g_bCl';

    % --- I_Ca ---
    parameter_names{12} = 'g_CaL';

    % --- I_NCX ---
    parameter_names{13} = 'Kdact';
    parameter_names{14} = 'KmCai';
    parameter_names{15} = 'KmCao';
    parameter_names{16} = 'KmNai';
    parameter_names{17} = 'KmNao';
    parameter_names{18} = 'g_NaCa';
    parameter_names{19} = 'ksat';
    parameter_names{20} = 'nu';

    % --- I_pCa ---
    parameter_names{21} = 'KmPCa';
    parameter_names{22} = 'g_pCa';

    % --- I_bCa ---
    parameter_names{23} = 'g_bCa';

    % --- I_f ---
    parameter_names{24} = 'E_f';
    parameter_names{25} = 'g_f';

    % --- Na Concentrations ---
    parameter_names{26} = 'Nao';

    % --- K Concentration ---
    parameter_names{27} = 'K_i';
    parameter_names{28} = 'Ko';

    % --- Ca Concentrations ---
    parameter_names{29} = 'ce';

    % --- RyRs ---
    parameter_names{30} = 'K_RyR';
    parameter_names{31} = 'alpha_RyR';
    parameter_names{32} = 'beta_RyR';
    parameter_names{33} = 'eta_RyR';
    parameter_names{34} = 'gamma_RyR';
    parameter_names{35} = 'lambda_RyR';

    % --- Intracellular volumes ---
    parameter_names{36} = 'Vc';
    parameter_names{37} = 'Vd';
    parameter_names{38} = 'Vn';
    parameter_names{39} = 'Vs';
    parameter_names{40} = 'Vsl';

    % --- SERCA pump ---
    parameter_names{41} = 'J_SERCA_bar';
    parameter_names{42} = 'K_c';
    parameter_names{43} = 'K_n';

    % --- Ca Buffers ---
    parameter_names{44} = 'B_tot_c';
    parameter_names{45} = 'B_tot_d';
    parameter_names{46} = 'B_tot_s';
    parameter_names{47} = 'B_tot_sl';
    parameter_names{48} = 'k_off_c';
    parameter_names{49} = 'k_off_d';
    parameter_names{50} = 'k_off_s';
    parameter_names{51} = 'k_off_sl';
    parameter_names{52} = 'k_on_c';
    parameter_names{53} = 'k_on_d';
    parameter_names{54} = 'k_on_s';
    parameter_names{55} = 'k_on_sl';

    % --- Ca Fluxes ---
    parameter_names{56} = 'alpha_d_c';
    parameter_names{57} = 'alpha_n_s';
    parameter_names{58} = 'alpha_sl_c';

    % --- Cl Concentrations ---
    parameter_names{59} = 'Cli';
    parameter_names{60} = 'Clo';

    % --- Membrane potential ---
    parameter_names{61} = 'Cm';
    parameter_names{62} = 'Frdy';
    parameter_names{63} = 'R';
    parameter_names{64} = 'Temp';
    parameter_names{65} = 'chi';
    parameter_names{66} = 'stim_amplitude';
    parameter_names{67} = 'stim_duration';
    parameter_names{68} = 'stim_period';
    parameter_names{69} = 'stim_start';

    % --- Mechanics ---
    parameter_names{70} = 'Fa';
    parameter_names{71} = 'H_mech';
    parameter_names{72} = 'K_mech';
    parameter_names{73} = 'Trop_frac';
    parameter_names{74} = 'a_mech';
    parameter_names{75} = 'b_mech';
    parameter_names{76} = 'kDX0';
    parameter_names{77} = 'kPD0';
    parameter_names{78} = 'kPX0';
    parameter_names{79} = 'kRD0';
    parameter_names{80} = 'kXD0';
    parameter_names{81} = 'kXP0';
    parameter_names{82} = 'k_off_t';
    parameter_names{83} = 'k_on_t';
    parameter_names{84} = 'kv_mech';
    parameter_names{85} = 'l0';
    parameter_names{86} = 'm_mech';
    varargout(1) = {parameter_names};
  end
end
