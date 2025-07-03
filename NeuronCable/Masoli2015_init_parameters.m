function [parameters, varargout] = Masoli2015_init_parameters()
  % % Default parameter values for ODE model: Masoli2015
  % % --------------------------------------------------
  % %
  % % parameters = Masoli2015_init_parameters();
  % % [parameters, parameters_names] = Masoli2015_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(33, 1);

  % --- Membrane ---
  parameters(1) = 0.001; % Cm;
  parameters(2) = 96485.3365; % F;
  parameters(3) = 96.4853365; % FARADAY;
  parameters(4) = 8.31446261815; % R;
  parameters(5) = 0.0; % Tdiff;
  parameters(6) = 37; % Temp;
  parameters(7) = 4.5e-05; % cai;
  parameters(8) = 2.0; % cao;

  % --- Stimulation ---
  parameters(9) = 0; % stim_amplitude;
  parameters(10) = 0.1; % stim_duration;
  parameters(11) = 800000; % stim_end;
  parameters(12) = 120.0; % stim_frequency;
  parameters(13) = 0; % stim_start;

  % --- Nav1.6 ---
  parameters(14) = 60; % E_Na;
  parameters(15) = 0.214; % g_Na;

  % --- Kv1.1 ---
  parameters(16) = -88; % E_K;
  parameters(17) = 0.002; % g_Kv11;

  % --- Kv1.5 ---
  parameters(18) = 0; % g_Kv15;

  % --- Kv3.3 ---
  parameters(19) = 0; % g_Kv33;

  % --- Kv3.4 ---
  parameters(20) = 0.05; % g_Kv34;

  % --- Kv4.3 ---
  parameters(21) = 0; % g_Kv43;

  % --- Kir2.x ---
  parameters(22) = 3e-05; % g_Kir2x;

  % --- Kca1.1 ---
  parameters(23) = 0.01; % g_Kca11;

  % --- Kca2.2 ---
  parameters(24) = 0.001; % g_Kca22;

  % --- Kca3.1 ---
  parameters(25) = 0.01; % g_Kca31;

  % --- Cav2.1 ---
  parameters(26) = 0.00022; % g_Cav21;

  % --- Cav3.1 ---
  parameters(27) = 7e-06; % g_Cav31;

  % --- Cav3.2 ---
  parameters(28) = 0.0008; % g_Cav32;

  % --- Cav3.3 ---
  parameters(29) = 0.0001; % g_Cav33;

  % --- HCN1 ---
  parameters(30) = -34.4; % E_h;
  parameters(31) = 0.0004; % g_HCN1;

  % --- leak ---
  parameters(32) = -63; % E_leak;
  parameters(33) = 0.0011; % g_leak;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(33, 1);

    % --- Membrane ---
    parameter_names{1} = 'Cm';
    parameter_names{2} = 'F';
    parameter_names{3} = 'FARADAY';
    parameter_names{4} = 'R';
    parameter_names{5} = 'Tdiff';
    parameter_names{6} = 'Temp';
    parameter_names{7} = 'cai';
    parameter_names{8} = 'cao';

    % --- Stimulation ---
    parameter_names{9} = 'stim_amplitude';
    parameter_names{10} = 'stim_duration';
    parameter_names{11} = 'stim_end';
    parameter_names{12} = 'stim_frequency';
    parameter_names{13} = 'stim_start';

    % --- Nav1.6 ---
    parameter_names{14} = 'E_Na';
    parameter_names{15} = 'g_Na';

    % --- Kv1.1 ---
    parameter_names{16} = 'E_K';
    parameter_names{17} = 'g_Kv11';

    % --- Kv1.5 ---
    parameter_names{18} = 'g_Kv15';

    % --- Kv3.3 ---
    parameter_names{19} = 'g_Kv33';

    % --- Kv3.4 ---
    parameter_names{20} = 'g_Kv34';

    % --- Kv4.3 ---
    parameter_names{21} = 'g_Kv43';

    % --- Kir2.x ---
    parameter_names{22} = 'g_Kir2x';

    % --- Kca1.1 ---
    parameter_names{23} = 'g_Kca11';

    % --- Kca2.2 ---
    parameter_names{24} = 'g_Kca22';

    % --- Kca3.1 ---
    parameter_names{25} = 'g_Kca31';

    % --- Cav2.1 ---
    parameter_names{26} = 'g_Cav21';

    % --- Cav3.1 ---
    parameter_names{27} = 'g_Cav31';

    % --- Cav3.2 ---
    parameter_names{28} = 'g_Cav32';

    % --- Cav3.3 ---
    parameter_names{29} = 'g_Cav33';

    % --- HCN1 ---
    parameter_names{30} = 'E_h';
    parameter_names{31} = 'g_HCN1';

    % --- leak ---
    parameter_names{32} = 'E_leak';
    parameter_names{33} = 'g_leak';
    varargout(1) = {parameter_names};
  end
end