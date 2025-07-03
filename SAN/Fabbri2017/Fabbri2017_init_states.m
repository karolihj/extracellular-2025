function [states, varargout] = Fabbri2017_init_states()
  % % Default state values for ODE model:
  % Human_SAN_Fabbri_Fantini_Wilders_Severi_2017


  % --- Default initial state values --- 
  states = zeros(33, 1);

  % --- FCa gate ---
  states(1) = 0.86423078518581; % fCa;

  % --- Ca SR release ---
  states(2) = 0.92169128697230; % R;
  states(3) = 0.00000000391251; % O;
  states(4) = 0.00000000033238; % I;
  states(5) = 0.07830771519954; % RI;

  % --- Ca buffering ---
  states(6) = 0.01715888580237; % fTC;
  states(7) = 0.25310992753848; % fTMC;
  states(8) = 0.65981618340199; % fTMM;
  states(9) = 0.20978972069401; % fCMi;
  states(10) = 0.13862767118132; % fCMs;
  states(11) = 0.12198138991798; % fCQ;

  % --- y gate ---
  states(12) = 0.01154178377492; % y;

  % --- m gate ---
  states(13) = 0.30499257534567; % m;

  % --- h gate ---
  states(14) = 0.00953685987234; % h;

  % --- DL gate ---
  states(15) = 0.00056014499225; % dL;

  % --- FL gate ---
  states(16) = 0.91694879510376; % fL;

  % --- DT gate ---
  states(17) = 0.12707870505249; % dT;

  % --- FT gate ---
  states(18) = 0.07337071841528; % fT;

  % --- Ca dynamics ---
  states(19) = 0.00008753086716; % Cai;
  states(20) = 0.00005309725840; % Ca_sub;
  states(21) = 0.37585331454361; % Ca_nsr;
  states(22) = 0.35265242438650; % Ca_jsr;

  % --- RKur gate ---
  states(23) = 0.00665591808436; % r_Kur;

  % --- SKur gate ---
  states(24) = 0.87447372646596; % s_Kur;

  % --- q gate ---
  states(25) = 0.51528094534181; % q;

  % --- r gate ---
  states(26) = 0.01044346002321; % r;

  % --- Pa gate ---
  states(27) = 0.28033464673098; % paS;
  states(28) = 0.00599041194084; % paF;

  % --- Pi gate ---
  states(29) = 0.76640584882757; % piy;

  % --- n gate ---
  states(30) = 0.09782669210473; % n;

  % --- a gate ---
  states(31) = 0.00290118067236; % a;

  % --- Nai_concentration ---
  states(32) = 5.00000000000000; % Nai_;

  % --- Membrane ---
  states(33) = -48.89158069082388; % V_ode;

  if nargout == 2

    % --- State names --- 
    state_names = cell(33, 1);

    % --- FCa gate ---
    state_names{1} = 'fCa';

    % --- Ca SR release ---
    state_names{2} = 'R';
    state_names{3} = 'O';
    state_names{4} = 'I';
    state_names{5} = 'RI';

    % --- Ca buffering ---
    state_names{6} = 'fTC';
    state_names{7} = 'fTMC';
    state_names{8} = 'fTMM';
    state_names{9} = 'fCMi';
    state_names{10} = 'fCMs';
    state_names{11} = 'fCQ';

    % --- y gate ---
    state_names{12} = 'y';

    % --- m gate ---
    state_names{13} = 'm';

    % --- h gate ---
    state_names{14} = 'h';

    % --- DL gate ---
    state_names{15} = 'dL';

    % --- FL gate ---
    state_names{16} = 'fL';

    % --- DT gate ---
    state_names{17} = 'dT';

    % --- FT gate ---
    state_names{18} = 'fT';

    % --- Ca dynamics ---
    state_names{19} = 'Cai';
    state_names{20} = 'Ca_sub';
    state_names{21} = 'Ca_nsr';
    state_names{22} = 'Ca_jsr';

    % --- RKur gate ---
    state_names{23} = 'r_Kur';

    % --- SKur gate ---
    state_names{24} = 's_Kur';

    % --- q gate ---
    state_names{25} = 'q';

    % --- r gate ---
    state_names{26} = 'r';

    % --- Pa gate ---
    state_names{27} = 'paS';
    state_names{28} = 'paF';

    % --- Pi gate ---
    state_names{29} = 'piy';

    % --- n gate ---
    state_names{30} = 'n';

    % --- a gate ---
    state_names{31} = 'a';

    % --- Nai_concentration ---
    state_names{32} = 'Nai_';

    % --- Membrane ---
    state_names{33} = 'V_ode';
    varargout(1) = {state_names};
  end
end