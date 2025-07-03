function [states, varargout] = Riz2014_init_states_metabolic()
  % % Default state values for ODE model: Riz2014 metabiloc version
  % % -------------------------------------------
  % %
  % % states = Riz2014_init_states_metabolic();
  % % [states, states_names] = Riz2014_init_states_metabolic();

  % --- Default initial state values --- 
  states = zeros(14, 1);

  % --- I_Na ---
  states(1) = 0.83; % h_Na;

  % --- I_CaL ---
  states(2) = 0.98; % h_CaL;

  % --- I_CaT ---
  states(3) = 0.17; % h_CaT;

  % --- I_BK ---
  states(4) = 0.0058; % m_BK;

  % --- I_Kv ---
  states(5) = 0.0056; % m_Kv;

  % --- I_HERG ---
  states(6) = 0.098; % m_HERG;
  states(7) = 0.64; % h_HERG;

  % --- Metabolism ---
  states(8) = 3.5; % G6PF6P;
  states(9) = 0.44; % FBP;
  states(10) = 0.028; % DHAPG3P;
  states(11) = 1.05; % a;

  % --- Calcium concentrations ---
  states(12) = 0.12; % Ca_m;
  states(13) = 0.067; % Ca_c;

  % --- Membrane potential ---
  states(14) = -51.4; % V;

  if nargout == 2

    % --- State names --- 
    state_names = cell(14, 1);

    % --- I_Na ---
    state_names{1} = 'h_Na';

    % --- I_CaL ---
    state_names{2} = 'h_CaL';

    % --- I_CaT ---
    state_names{3} = 'h_CaT';

    % --- I_BK ---
    state_names{4} = 'm_BK';

    % --- I_Kv ---
    state_names{5} = 'm_Kv';

    % --- I_HERG ---
    state_names{6} = 'm_HERG';
    state_names{7} = 'h_HERG';

    % --- Metabolism ---
    state_names{8} = 'G6PF6P';
    state_names{9} = 'FBP';
    state_names{10} = 'DHAPG3P';
    state_names{11} = 'a';

    % --- Calcium concentrations ---
    state_names{12} = 'Ca_m';
    state_names{13} = 'Ca_c';

    % --- Membrane potential ---
    state_names{14} = 'V';
    varargout(1) = {state_names};
  end
end
