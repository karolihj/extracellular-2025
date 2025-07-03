function [parameters, varargout] = Fabbri2017_init_parameters()
  % % Default parameter values for ODE model:
  % Human_SAN_Fabbri_Fantini_Wilders_Severi_2017

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(91, 1);

  % --- Membrane ---
  parameters(1) = 5.7e-05; % C;
  parameters(2) = 96485.3415; % F;
  parameters(3) = 8314.472; % R_Membrane;
  parameters(4) = 310.0; % T;
  parameters(5) = 0; % clamp_mode;

  % --- Voltage clamp ---
  parameters(6) = -45.0; % V_holding;
  parameters(7) = -35.0; % V_test;
  parameters(8) = 0.5; % t_holding;
  parameters(9) = 0.5; % t_test;

  % --- Rate modulation experiments ---
  parameters(10) = 0.0; % ACh;
  parameters(11) = 0; % Iso_1_uM;

  % --- i_CaL ---
  parameters(12) = 0.4578; % P_CaL;

  % --- FCa gate ---
  parameters(13) = 0.000338; % Km_fCa;
  parameters(14) = 0.0075; % alpha_fCa;

  % --- DL gate ---
  parameters(15) = -16.4508; % V_dL;
  parameters(16) = 4.3371; % k_dL;

  % --- FL gate ---
  parameters(17) = 0.0; % k_fL;
  parameters(18) = 0.0; % shift_fL;

  % --- Ca SR release ---
  parameters(19) = 0.45; % EC50_SR;
  parameters(20) = 2.5; % HSR;
  parameters(21) = 15; % MaxSR;
  parameters(22) = 1; % MinSR;
  parameters(23) = 500.0; % kiCa;
  parameters(24) = 5.0; % kim;
  parameters(25) = 10000.0; % koCa;
  parameters(26) = 660.0; % kom;
  parameters(27) = 148041085.1; % ks;

  % --- Ca intracellular fluxes ---
  parameters(28) = 0.000286113; % K_up;
  parameters(29) = 5.0; % P_up_basal;
  parameters(30) = 5e-05; % slope_up;
  parameters(31) = 5.469e-05; % tau_dif_Ca;
  parameters(32) = 0.04; % tau_tr;

  % --- Ca buffering ---
  parameters(33) = 0.045; % CM_tot;
  parameters(34) = 10.0; % CQ_tot;
  parameters(35) = 2.5; % Mgi;
  parameters(36) = 0.031; % TC_tot;
  parameters(37) = 0.062; % TMC_tot;
  parameters(38) = 542.0; % kb_CM;
  parameters(39) = 445.0; % kb_CQ;
  parameters(40) = 446.0; % kb_TC;
  parameters(41) = 7.51; % kb_TMC;
  parameters(42) = 751.0; % kb_TMM;
  parameters(43) = 1642000.0; % kf_CM;
  parameters(44) = 175.4; % kf_CQ;
  parameters(45) = 88800.0; % kf_TC;
  parameters(46) = 227700.0; % kf_TMC;
  parameters(47) = 2277.0; % kf_TMM;

  % --- Cell parameters ---
  parameters(48) = 67.0; % L_cell;
  parameters(49) = 0.02; % L_sub;
  parameters(50) = 3.9; % R_cell;
  parameters(51) = 0.46; % V_i_part;
  parameters(52) = 0.0012; % V_jsr_part;
  parameters(53) = 0.0116; % V_nsr_part;

  % --- Ionic values ---
  parameters(54) = 1.8; % Cao;
  parameters(55) = 140.0; % Ki;
  parameters(56) = 5.4; % Ko;
  parameters(57) = 140.0; % Nao;

  % --- Nai_concentration ---
  parameters(58) = 1; % Nai_clamp;

  % --- i_f ---
  parameters(59) = 45.0; % Km_f;
  parameters(60) = 0.5927; % alpha;
  parameters(61) = 0; % blockade;
  parameters(62) = 0.00427; % g_f;

  % --- y gate ---
  parameters(63) = 0.0; % y_shift;

  % --- i_NaK ---
  parameters(64) = 1.4; % Km_Kp;
  parameters(65) = 14.0; % Km_Nap;
  parameters(66) = 0.08105; % i_NaK_max;

  % --- i_NaCa ---
  parameters(67) = 395.3; % K1ni;
  parameters(68) = 1628.0; % K1no;
  parameters(69) = 2.289; % K2ni;
  parameters(70) = 561.4; % K2no;
  parameters(71) = 26.44; % K3ni;
  parameters(72) = 4.663; % K3no;
  parameters(73) = 3.343; % K_NaCa;
  parameters(74) = 0.0207; % Kci;
  parameters(75) = 26.44; % Kcni;
  parameters(76) = 3.663; % Kco;
  parameters(77) = 0.1369; % Qci;
  parameters(78) = 0; % Qco;
  parameters(79) = 0.4315; % Qn;
  parameters(80) = 0; % blockade_NaCa;

  % --- i_Na ---
  parameters(81) = 0.0223; % g_Na;
  parameters(82) = 0.0; % g_Na_L;

  % --- m gate ---
  parameters(83) = 1e-05; % delta_m;

  % --- i_CaT ---
  parameters(84) = 0.04132; % P_CaT;

  % --- FT gate ---
  parameters(85) = 0.0; % offset_fT;

  % --- i_Kur ---
  parameters(86) = 0.0001539; % g_Kur;

  % --- i_to ---
  parameters(87) = 0.0035; % g_to;

  % --- i_Kr ---
  parameters(88) = 0.00424; % g_Kr;

  % --- i_Ks ---
  parameters(89) = 0.00065; % g_Ks_;

  % --- i_KACh ---
  parameters(90) = 1; % ACh_on;
  parameters(91) = 0.00345; % g_KACh;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(91, 1);

    % --- Membrane ---
    parameter_names{1} = 'C';
    parameter_names{2} = 'F';
    parameter_names{3} = 'R_Membrane';
    parameter_names{4} = 'T';
    parameter_names{5} = 'clamp_mode';

    % --- Voltage clamp ---
    parameter_names{6} = 'V_holding';
    parameter_names{7} = 'V_test';
    parameter_names{8} = 't_holding';
    parameter_names{9} = 't_test';

    % --- Rate modulation experiments ---
    parameter_names{10} = 'ACh';
    parameter_names{11} = 'Iso_1_uM';

    % --- i_CaL ---
    parameter_names{12} = 'P_CaL';

    % --- FCa gate ---
    parameter_names{13} = 'Km_fCa';
    parameter_names{14} = 'alpha_fCa';

    % --- DL gate ---
    parameter_names{15} = 'V_dL';
    parameter_names{16} = 'k_dL';

    % --- FL gate ---
    parameter_names{17} = 'k_fL';
    parameter_names{18} = 'shift_fL';

    % --- Ca SR release ---
    parameter_names{19} = 'EC50_SR';
    parameter_names{20} = 'HSR';
    parameter_names{21} = 'MaxSR';
    parameter_names{22} = 'MinSR';
    parameter_names{23} = 'kiCa';
    parameter_names{24} = 'kim';
    parameter_names{25} = 'koCa';
    parameter_names{26} = 'kom';
    parameter_names{27} = 'ks';

    % --- Ca intracellular fluxes ---
    parameter_names{28} = 'K_up';
    parameter_names{29} = 'P_up_basal';
    parameter_names{30} = 'slope_up';
    parameter_names{31} = 'tau_dif_Ca';
    parameter_names{32} = 'tau_tr';

    % --- Ca buffering ---
    parameter_names{33} = 'CM_tot';
    parameter_names{34} = 'CQ_tot';
    parameter_names{35} = 'Mgi';
    parameter_names{36} = 'TC_tot';
    parameter_names{37} = 'TMC_tot';
    parameter_names{38} = 'kb_CM';
    parameter_names{39} = 'kb_CQ';
    parameter_names{40} = 'kb_TC';
    parameter_names{41} = 'kb_TMC';
    parameter_names{42} = 'kb_TMM';
    parameter_names{43} = 'kf_CM';
    parameter_names{44} = 'kf_CQ';
    parameter_names{45} = 'kf_TC';
    parameter_names{46} = 'kf_TMC';
    parameter_names{47} = 'kf_TMM';

    % --- Cell parameters ---
    parameter_names{48} = 'L_cell';
    parameter_names{49} = 'L_sub';
    parameter_names{50} = 'R_cell';
    parameter_names{51} = 'V_i_part';
    parameter_names{52} = 'V_jsr_part';
    parameter_names{53} = 'V_nsr_part';

    % --- Ionic values ---
    parameter_names{54} = 'Cao';
    parameter_names{55} = 'Ki';
    parameter_names{56} = 'Ko';
    parameter_names{57} = 'Nao';

    % --- Nai_concentration ---
    parameter_names{58} = 'Nai_clamp';

    % --- i_f ---
    parameter_names{59} = 'Km_f';
    parameter_names{60} = 'alpha';
    parameter_names{61} = 'blockade';
    parameter_names{62} = 'g_f';

    % --- y gate ---
    parameter_names{63} = 'y_shift';

    % --- i_NaK ---
    parameter_names{64} = 'Km_Kp';
    parameter_names{65} = 'Km_Nap';
    parameter_names{66} = 'i_NaK_max';

    % --- i_NaCa ---
    parameter_names{67} = 'K1ni';
    parameter_names{68} = 'K1no';
    parameter_names{69} = 'K2ni';
    parameter_names{70} = 'K2no';
    parameter_names{71} = 'K3ni';
    parameter_names{72} = 'K3no';
    parameter_names{73} = 'K_NaCa';
    parameter_names{74} = 'Kci';
    parameter_names{75} = 'Kcni';
    parameter_names{76} = 'Kco';
    parameter_names{77} = 'Qci';
    parameter_names{78} = 'Qco';
    parameter_names{79} = 'Qn';
    parameter_names{80} = 'blockade_NaCa';

    % --- i_Na ---
    parameter_names{81} = 'g_Na';
    parameter_names{82} = 'g_Na_L';

    % --- m gate ---
    parameter_names{83} = 'delta_m';

    % --- i_CaT ---
    parameter_names{84} = 'P_CaT';

    % --- FT gate ---
    parameter_names{85} = 'offset_fT';

    % --- i_Kur ---
    parameter_names{86} = 'g_Kur';

    % --- i_to ---
    parameter_names{87} = 'g_to';

    % --- i_Kr ---
    parameter_names{88} = 'g_Kr';

    % --- i_Ks ---
    parameter_names{89} = 'g_Ks_';

    % --- i_KACh ---
    parameter_names{90} = 'ACh_on';
    parameter_names{91} = 'g_KACh';
    varargout(1) = {parameter_names};
  end
end