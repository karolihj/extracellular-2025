function [states, varargout] = Masoli2015_init_states()
  % % Default state values for ODE model: Masoli2015
  % % ----------------------------------------------
  % %
  % % states = Masoli2015_init_states();
  % % [states, states_names] = Masoli2015_init_states();

  % --- Default initial state values --- 
  states = zeros(46, 1);

  % --- Nav1.6 ---
  states(1) = 0.745657349519; % C1_Na;
  states(2) = 0.18681749783; % C2_Na;
  states(3) = 0.0175508753271; % C3_Na;
  states(4) = 0.00073259055439; % C4_Na;
  states(5) = 1.136427441e-05; % C5_Na;
  states(6) = 4.176150917e-05; % O_Na;
  states(7) = 0.00012203268796; % B_Na;
  states(8) = 0.0054257348267; % I1_Na;
  states(9) = 0.0150081485653; % I2_Na;
  states(10) = 0.0155777872341; % I3_Na;
  states(11) = 0.00718538355943; % I4_Na;
  states(12) = 0.00123878700626; % I5_Na;

  % --- Kv1.1 ---
  states(13) = 0.0784327349458; % n_Kv11;

  % --- Kv1.5 ---
  states(14) = 0.062693414843; % m_Kv15;
  states(15) = 0.953843893364; % n_Kv15;
  states(16) = 0.982440214681; % u_Kv15;

  % --- Kv3.3 ---
  states(17) = 0.0206633026411; % n_Kv33;

  % --- Kv3.4 ---
  states(18) = 0.143344433759; % m_Kv34;
  states(19) = 0.882087088068; % h_Kv34;

  % --- Kv4.3 ---
  states(20) = 0.0299743553258; % a_Kv43;
  states(21) = 0.185952387828; % b_Kv43;

  % --- Kir2.x ---
  states(22) = 0.206778862228; % d_Kir2x;

  % --- Kca1.1 ---
  states(23) = 0.983767260817; % C0_Kca11;
  states(24) = 0.0161123619009; % C1_Kca11;
  states(25) = 9.895934853e-05; % C2_Kca11;
  states(26) = 2.7012724e-07; % C3_Kca11;
  states(27) = 1.801623e-05; % O0_Kca11;
  states(28) = 2.94586477e-06; % O1_Kca11;
  states(29) = 1.8046753e-07; % O2_Kca11;
  states(30) = 4.91765e-09; % O3_Kca11;
  states(31) = 5.012e-11; % O4_Kca11;

  % --- Kca2.2 ---
  states(32) = 0.88748612049; % c1_Kca22;
  states(33) = 0.0999548330279; % c2_Kca22;
  states(34) = 0.00900749464037; % c3_Kca22;
  states(35) = 0.00144123766793; % o1_Kca22;
  states(36) = 0.00194798850465; % o2_Kca22;

  % --- Kca3.1 ---
  states(37) = 0.00342250393833; % Y_Kca31;

  % --- Cav2.1 ---
  states(38) = 0.0113511215154; % m_Cav21;

  % --- Cav3.1 ---
  states(39) = 0.176172416128; % m_Cav31;
  states(40) = 0.191716948331; % h_Cav31;

  % --- Cav3.2 ---
  states(41) = 0.15549099887; % m_Cav32;
  states(42) = 0.0574480599015; % h_Cav32;

  % --- Cav3.3 ---
  states(43) = 0.440000580693; % n_Cav33;
  states(44) = 0.164229960931; % l_Cav33;

  % --- HCN1 ---
  states(45) = 0.00144283978473; % h_HCN1;

  % --- Membrane ---
  states(46) = -72.9572879674; % v;

  if nargout == 2

    % --- State names --- 
    state_names = cell(46, 1);

    % --- Nav1.6 ---
    state_names{1} = 'C1_Na';
    state_names{2} = 'C2_Na';
    state_names{3} = 'C3_Na';
    state_names{4} = 'C4_Na';
    state_names{5} = 'C5_Na';
    state_names{6} = 'O_Na';
    state_names{7} = 'B_Na';
    state_names{8} = 'I1_Na';
    state_names{9} = 'I2_Na';
    state_names{10} = 'I3_Na';
    state_names{11} = 'I4_Na';
    state_names{12} = 'I5_Na';

    % --- Kv1.1 ---
    state_names{13} = 'n_Kv11';

    % --- Kv1.5 ---
    state_names{14} = 'm_Kv15';
    state_names{15} = 'n_Kv15';
    state_names{16} = 'u_Kv15';

    % --- Kv3.3 ---
    state_names{17} = 'n_Kv33';

    % --- Kv3.4 ---
    state_names{18} = 'm_Kv34';
    state_names{19} = 'h_Kv34';

    % --- Kv4.3 ---
    state_names{20} = 'a_Kv43';
    state_names{21} = 'b_Kv43';

    % --- Kir2.x ---
    state_names{22} = 'd_Kir2x';

    % --- Kca1.1 ---
    state_names{23} = 'C0_Kca11';
    state_names{24} = 'C1_Kca11';
    state_names{25} = 'C2_Kca11';
    state_names{26} = 'C3_Kca11';
    state_names{27} = 'O0_Kca11';
    state_names{28} = 'O1_Kca11';
    state_names{29} = 'O2_Kca11';
    state_names{30} = 'O3_Kca11';
    state_names{31} = 'O4_Kca11';

    % --- Kca2.2 ---
    state_names{32} = 'c1_Kca22';
    state_names{33} = 'c2_Kca22';
    state_names{34} = 'c3_Kca22';
    state_names{35} = 'o1_Kca22';
    state_names{36} = 'o2_Kca22';

    % --- Kca3.1 ---
    state_names{37} = 'Y_Kca31';

    % --- Cav2.1 ---
    state_names{38} = 'm_Cav21';

    % --- Cav3.1 ---
    state_names{39} = 'm_Cav31';
    state_names{40} = 'h_Cav31';

    % --- Cav3.2 ---
    state_names{41} = 'm_Cav32';
    state_names{42} = 'h_Cav32';

    % --- Cav3.3 ---
    state_names{43} = 'n_Cav33';
    state_names{44} = 'l_Cav33';

    % --- HCN1 ---
    state_names{45} = 'h_HCN1';

    % --- Membrane ---
    state_names{46} = 'v';
    varargout(1) = {state_names};
  end
end