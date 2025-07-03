function [states, varargout] = base_model_init_states(filename)
  % % Default state values for ODE model: base_model

  load(filename, 'states')

  if nargout == 2

    % --- State names --- 
    state_names = cell(32, 1);

    % --- I_Na ---
    state_names{1} = 'm';
    state_names{2} = 'j';

    % --- I_NaL ---
    state_names{3} = 'mL';
    state_names{4} = 'hL';

    % --- I_Kr ---
    state_names{5} = 'Xr1';
    state_names{6} = 'Xr2';

    % --- I_Ks ---
    state_names{7} = 'x_Ks';

    % --- i_to ---
    state_names{8} = 'q';
    state_names{9} = 'r';

    % --- I_Ca ---
    state_names{10} = 'd';
    state_names{11} = 'f';
    state_names{12} = 'f_Ca_B';

    % --- I_f ---
    state_names{13} = 'xf';

    % --- RyRs ---
    state_names{14} = 'r_RyR';

    % --- Mechanics ---
    state_names{15} = 'bt';
    state_names{16} = 'R_mech';
    state_names{17} = 'D_mech';
    state_names{18} = 'X_mech';
    state_names{19} = 'P_mech';
    state_names{20} = 'l';
    state_names{21} = 'lv';

    % --- Ca Concentrations ---
    state_names{22} = 'cn';
    state_names{23} = 'cc';
    state_names{24} = 'cd';
    state_names{25} = 'csl';
    state_names{26} = 'cs';

    % --- Ca Buffer Concentrations ---
    state_names{27} = 'bc';
    state_names{28} = 'bd';
    state_names{29} = 'bs';
    state_names{30} = 'bsl';

    % --- Membrane potential ---
    state_names{31} = 'V_m';

    % --- Na Concentrations ---
    state_names{32} = 'Na_i';
    varargout(1) = {state_names};
  end
end
