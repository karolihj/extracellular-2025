function [values] = system_rhs(t, states, P)
  % Compute the right hand side of the Masoli2015 ODE

  N = P.N;
  parameters = P.parameters;
  G = P.G;

  Ns = 46;
  Np = 33;


  % Assign states
  C1_Na=states(1:Ns:Ns*N); C2_Na=states(2:Ns:Ns*N); C3_Na=states(3:Ns:Ns*N); C4_Na=states(4:Ns:Ns*N);...
    C5_Na=states(5:Ns:Ns*N); O_Na=states(6:Ns:Ns*N); B_Na=states(7:Ns:Ns*N); I1_Na=states(8:Ns:Ns*N);...
    I2_Na=states(9:Ns:Ns*N); I3_Na=states(10:Ns:Ns*N); I4_Na=states(11:Ns:Ns*N); I5_Na=states(12:Ns:Ns*N);...
    n_Kv11=states(13:Ns:Ns*N); m_Kv15=states(14:Ns:Ns*N); n_Kv15=states(15:Ns:Ns*N);...
    u_Kv15=states(16:Ns:Ns*N); n_Kv33=states(17:Ns:Ns*N); m_Kv34=states(18:Ns:Ns*N);...
    h_Kv34=states(19:Ns:Ns*N); a_Kv43=states(20:Ns:Ns*N); b_Kv43=states(21:Ns:Ns*N);...
    d_Kir2x=states(22:Ns:Ns*N); C0_Kca11=states(23:Ns:Ns*N); C1_Kca11=states(24:Ns:Ns*N);...
    C2_Kca11=states(25:Ns:Ns*N); C3_Kca11=states(26:Ns:Ns*N); O0_Kca11=states(27:Ns:Ns*N);...
    O1_Kca11=states(28:Ns:Ns*N); O2_Kca11=states(29:Ns:Ns*N); O3_Kca11=states(30:Ns:Ns*N);...
    O4_Kca11=states(31:Ns:Ns*N); c1_Kca22=states(32:Ns:Ns*N); c2_Kca22=states(33:Ns:Ns*N);...
    c3_Kca22=states(34:Ns:Ns*N); o1_Kca22=states(35:Ns:Ns*N); o2_Kca22=states(36:Ns:Ns*N);...
    Y_Kca31=states(37:Ns:Ns*N); m_Cav21=states(38:Ns:Ns*N); m_Cav31=states(39:Ns:Ns*N);...
    h_Cav31=states(40:Ns:Ns*N); m_Cav32=states(41:Ns:Ns*N); h_Cav32=states(42:Ns:Ns*N);...
    n_Cav33=states(43:Ns:Ns*N); l_Cav33=states(44:Ns:Ns*N); h_HCN1=states(45:Ns:Ns*N); v=states(46:Ns:Ns*N);

  % Assign parameters
  Cm=parameters(1:Np:end); F=parameters(2:Np:end); FARADAY=parameters(3:Np:end); R=parameters(4:Np:end);...
    Tdiff=parameters(5:Np:end); Temp=parameters(6:Np:end); cai=parameters(7:Np:end);...
    cao=parameters(8:Np:end); stim_amplitude=parameters(9:Np:end);...
    stim_duration=parameters(10:Np:end); stim_frequency=parameters(12:Np:end);...
    stim_start=parameters(13:Np:end); E_Na=parameters(14:Np:end); g_Na=parameters(15:Np:end);...
    E_K=parameters(16:Np:end); g_Kv11=parameters(17:Np:end); g_Kv15=parameters(18:Np:end);...
    g_Kv33=parameters(19:Np:end); g_Kv34=parameters(20:Np:end); g_Kv43=parameters(21:Np:end);...
    g_Kir2x=parameters(22:Np:end); g_Kca11=parameters(23:Np:end); g_Kca22=parameters(24:Np:end);...
    g_Kca31=parameters(25:Np:end); g_Cav21=parameters(26:Np:end); g_Cav31=parameters(27:Np:end);...
    g_Cav32=parameters(28:Np:end); g_Cav33=parameters(29:Np:end); E_h=parameters(30:Np:end);...
    g_HCN1=parameters(31:Np:end); E_leak=parameters(32:Np:end); g_leak=parameters(33:Np:end);

  % Init return args
  values = zeros(size(states));

  % Expressions for the Nav1.6 component
  allow_change = ((t >= Tdiff).*(1) + ~(t >= Tdiff).*(0));
  q10_Na = 3;
  qt_Na = q10_Na.^(-2.2 + 0.1.*Temp);
  Con_Na = 0.005;
  Coff_Na = 0.5;
  Oon_Na = 0.75;
  Ooff_Na = 0.005;
  alpha_Na = 150;
  beta_Na = 3;
  gamma_Na = 150;
  delta_Na = 40;
  epsilon_Na = 1.75;
  zeta_Na = 0.03;
  x1_Na = 20;
  x2_Na = -20;
  x3_Na = 1e+12;
  x4_Na = -1e+12;
  x5_Na = 1e+12;
  x6_Na = -25;
  alfac_Na = (Oon_Na./Con_Na).^0.25;
  btfac_Na = (Ooff_Na./Coff_Na).^0.25;
  f01_Na = 4.*alpha_Na.*exp(v./x1_Na).*qt_Na;
  f02_Na = 3.*alpha_Na.*exp(v./x1_Na).*qt_Na;
  f03_Na = 2.*alpha_Na.*exp(v./x1_Na).*qt_Na;
  f04_Na = alpha_Na.*exp(v./x1_Na).*qt_Na;
  f0O_Na = gamma_Na.*exp(v./x3_Na).*qt_Na;
  fip_Na = epsilon_Na.*exp(v./x5_Na).*qt_Na;
  f11_Na = 4.*alpha_Na.*alfac_Na.*exp(v./x1_Na).*qt_Na;
  f12_Na = 3.*alpha_Na.*alfac_Na.*exp(v./x1_Na).*qt_Na;
  f13_Na = 2.*alpha_Na.*alfac_Na.*exp(v./x1_Na).*qt_Na;
  f14_Na = alpha_Na.*alfac_Na.*exp(v./x1_Na).*qt_Na;
  f1n_Na = gamma_Na.*exp(v./x3_Na).*qt_Na;
  fi1_Na = Con_Na.*qt_Na;
  fi2_Na = Con_Na.*alfac_Na.*qt_Na;
  fi3_Na = Con_Na.*alfac_Na.^2.*qt_Na;
  fi4_Na = Con_Na.*alfac_Na.^3.*qt_Na;
  fi5_Na = Con_Na.*alfac_Na.^4.*qt_Na;
  fin_Na = Oon_Na.*qt_Na;
  b01_Na = beta_Na.*exp(v./x2_Na).*qt_Na;
  b02_Na = 2.*beta_Na.*exp(v./x2_Na).*qt_Na;
  b03_Na = 3.*beta_Na.*exp(v./x2_Na).*qt_Na;
  b04_Na = 4.*beta_Na.*exp(v./x2_Na).*qt_Na;
  b0O_Na = delta_Na.*exp(v./x4_Na).*qt_Na;
  bip_Na = zeta_Na.*exp(v./x6_Na).*qt_Na;
  b11_Na = beta_Na.*btfac_Na.*exp(v./x2_Na).*qt_Na;
  b12_Na = 2.*beta_Na.*btfac_Na.*exp(v./x2_Na).*qt_Na;
  b13_Na = 3.*beta_Na.*btfac_Na.*exp(v./x2_Na).*qt_Na;
  b14_Na = 4.*beta_Na.*btfac_Na.*exp(v./x2_Na).*qt_Na;
  b1n_Na = delta_Na.*exp(v./x4_Na).*qt_Na;
  bi1_Na = Coff_Na.*qt_Na;
  bi2_Na = Coff_Na.*btfac_Na.*qt_Na;
  bi3_Na = Coff_Na.*btfac_Na.^2.*qt_Na;
  bi4_Na = Coff_Na.*btfac_Na.^3.*qt_Na;
  bi5_Na = Coff_Na.*btfac_Na.^4.*qt_Na;
  bin_Na = Ooff_Na.*qt_Na;
  I6_Na = 1 - B_Na - C1_Na - C2_Na - C3_Na - C4_Na - C5_Na - I1_Na - I2_Na -...
    I3_Na - I4_Na - I5_Na - O_Na;
  values(1:Ns:end) = ((-f01_Na - fi1_Na).*C1_Na + C2_Na.*b01_Na +...
    I1_Na.*bi1_Na).*allow_change;
  values(2:Ns:end) = ((-b01_Na - f02_Na - fi2_Na).*C2_Na + C1_Na.*f01_Na +...
    C3_Na.*b02_Na + I2_Na.*bi2_Na).*allow_change;
  values(3:Ns:end) = ((-b02_Na - f03_Na - fi3_Na).*C3_Na + C2_Na.*f02_Na +...
    C4_Na.*b03_Na + I3_Na.*bi3_Na).*allow_change;
  values(4:Ns:end) = ((-b03_Na - f04_Na - fi4_Na).*C4_Na + C3_Na.*f03_Na +...
    C5_Na.*b04_Na + I4_Na.*bi4_Na).*allow_change;
  values(5:Ns:end) = ((-b04_Na - f0O_Na - fi5_Na).*C5_Na + C4_Na.*f04_Na +...
    I5_Na.*bi5_Na + O_Na.*b0O_Na).*allow_change;
  values(6:Ns:end) = ((-b0O_Na - fin_Na - fip_Na).*O_Na + B_Na.*bip_Na + C5_Na.*f0O_Na...
    + I6_Na.*bin_Na).*allow_change;
  values(7:Ns:end) = (O_Na.*fip_Na - B_Na.*bip_Na).*allow_change;
  values(8:Ns:end) = ((-bi1_Na - f11_Na).*I1_Na + C1_Na.*fi1_Na +...
    I2_Na.*b11_Na).*allow_change;
  values(9:Ns:end) = ((-b11_Na - bi2_Na - f12_Na).*I2_Na + C2_Na.*fi2_Na +...
    I1_Na.*f11_Na + I3_Na.*b12_Na).*allow_change;
  values(10:Ns:end) = ((-b12_Na - bi3_Na - f13_Na).*I3_Na + C3_Na.*fi3_Na +...
    I2_Na.*f12_Na + I4_Na.*b13_Na).*allow_change;
  values(11:Ns:end) = ((-b13_Na - bi4_Na - f14_Na).*I4_Na + C4_Na.*fi4_Na +...
    I3_Na.*f13_Na + I5_Na.*b14_Na).*allow_change;
  values(12:Ns:end) = ((-b14_Na - bi5_Na - f1n_Na).*I5_Na + C5_Na.*fi5_Na +...
    I4_Na.*f14_Na + I6_Na.*b1n_Na).*allow_change;
  I_Na = g_Na.*(-E_Na + v).*O_Na;

  % Expressions for the Kv1.1 component
  q10_Kv11 = 2.7;
  ca_Kv11 = 0.12889;
  cva_Kv11 = 45;
  cka_Kv11 = -33.90877;
  cb_Kv11 = 0.12889;
  cvb_Kv11 = 45;
  ckb_Kv11 = 12.42101;
  qt_Kv11 = q10_Kv11.^(-2.2 + 0.1.*Temp);
  alphan_Kv11 = ca_Kv11.*exp((-cva_Kv11 - v)./cka_Kv11);
  betan_Kv11 = cb_Kv11.*exp((-cvb_Kv11 - v)./ckb_Kv11);
  ninf_Kv11 = alphan_Kv11./(alphan_Kv11 + betan_Kv11);
  taun_Kv11 = 1./((alphan_Kv11 + betan_Kv11).*qt_Kv11);
  values(13:Ns:end) = (-n_Kv11 + ninf_Kv11).*allow_change./taun_Kv11;
  I_Kv11 = g_Kv11.*n_Kv11.^4.*(-E_K + v);

  % Expressions for the Kv1.5 component
  q10_Kv15 = 2.2.^(-3.7 + 0.1.*Temp);
  am_Kv15 = 0.65.*q10_Kv15./(1.66275285675.*exp(-0.0169491525424.*v) +...
    0.308365167897.*exp(-0.117647058824.*v));
  bm_Kv15 = 0.65.*q10_Kv15./(2.5 + 124.40338764.*exp(0.0588235294118.*v));
  mtau_Kv15 = 1./(3.*(am_Kv15 + bm_Kv15));
  minf_Kv15 = 1.0./(1.0 + 0.0425851362888.*exp(-0.104166666667.*v));
  an_Kv15 = 0.001.*q10_Kv15./(2.4 + 3.43809189356.*exp(-0.0128205128205.*v));
  bn_Kv15 = 2.75364493497e-08.*exp(0.0625.*v).*q10_Kv15;
  ntau_Kv15 = 0.333333333333./(an_Kv15 + bn_Kv15);
  ninf_Kv15 = 0.25 + 1.0./(1.35 + 1.6487212707.*exp(0.0714285714286.*v));
  uinf_Kv15 = 0.1 + 1.0./(1.1 + 1.6487212707.*exp(0.0714285714286.*v));
  utau_Kv15 = 6800.0;
  values(14:Ns:end) = (-m_Kv15 + minf_Kv15).*allow_change./mtau_Kv15;
  values(15:Ns:end) = (-n_Kv15 + ninf_Kv15).*allow_change./ntau_Kv15;
  values(16:Ns:end) = (-u_Kv15 + uinf_Kv15).*allow_change./utau_Kv15;
  I_Kv15 = g_Kv15.*m_Kv15.^3.*(0.1 + 1.0./(1.0 +...
    3.17036319489.*exp(-0.0769230769231.*v))).*(-E_K + v).*n_Kv15.*u_Kv15;

  % Expressions for the Kv3.3 component
  q10_Kv33 = 2.7;
  qt_Kv33 = q10_Kv33.^(-2.2 + 0.1.*Temp);
  ca_Kv33 = 0.22;
  cva_Kv33 = 16;
  cka_Kv33 = -26.5;
  cb_Kv33 = 0.22;
  cvb_Kv33 = 16;
  ckb_Kv33 = 26.5;
  alpha_Kv33 = ca_Kv33.*exp((-cva_Kv33 - v)./cka_Kv33).*qt_Kv33;
  beta_Kv33 = cb_Kv33.*exp((-cvb_Kv33 - v)./ckb_Kv33).*qt_Kv33;
  values(17:Ns:end) = ((1 - n_Kv33).*alpha_Kv33 - beta_Kv33.*n_Kv33).*allow_change;
  I_Kv33 = g_Kv33.*n_Kv33.^4.*(-E_K + v);

  % Expressions for the Kv3.4 component
  q10_Kv34 = 3;
  qt_Kv34 = q10_Kv34.^(-3.7 + 0.1.*Temp);
  mivh_Kv34 = -24;
  mik_Kv34 = 15.4;
  mty0_Kv34 = 0.00012851;
  mtvh1_Kv34 = 100.7;
  mtk1_Kv34 = 12.9;
  mtvh2_Kv34 = -56.0;
  mtk2_Kv34 = -23.1;
  hiy0_Kv34 = 0.31;
  hiA_Kv34 = 0.69;
  hivh_Kv34 = -5.802;
  hik_Kv34 = 11.2;
  minf_Kv34 = 1.0./(1.0 + exp((-11.0 + mivh_Kv34 - v)./mik_Kv34));
  mtau_Kv34 = 1000.*((11 + v < -35).*(0.000102675 +...
    0.0220402902836.*exp(0.0353481795688.*v)) + ~(11 + v < -35).*(mty0_Kv34 +...
    1.0./(exp((11.0 + mtvh1_Kv34 + v)./mtk1_Kv34) + exp((11.0 + mtvh2_Kv34 +...
    v)./mtk2_Kv34))))./qt_Kv34;
  hinf_Kv34 = hiy0_Kv34 + hiA_Kv34./(1 + exp((11.0 - hivh_Kv34 + v)./hik_Kv34));
  htau_Kv34 = 1000.*((11 + v > 0).*(0.0012 + 0.000487682413466.*exp(-0.141.*v)) +...
    ~(11 + v > 0).*(1.2202e-05 + 0.012.*exp(-(1.35685483871 +...
    0.0201612903226.*v).^2)))./qt_Kv34;
  values(18:Ns:end) = (-m_Kv34 + minf_Kv34).*allow_change./mtau_Kv34;
  values(19:Ns:end) = (-h_Kv34 + hinf_Kv34).*allow_change./htau_Kv34;
  I_Kv34 = g_Kv34.*m_Kv34.^3.*(-E_K + v).*h_Kv34;

  % Expressions for the Kv4.3 component
  Q10_Kv43 = 3.^(-2.55 + 0.1.*Temp);
  Aalpha_a_Kv43 = 0.8147;
  Kalpha_a_Kv43 = -23.32708;
  V0alpha_a_Kv43 = -9.17203;
  Abeta_a_Kv43 = 0.1655;
  Kbeta_a_Kv43 = 19.47175;
  V0beta_a_Kv43 = -18.27914;
  Aalpha_b_Kv43 = 0.0368;
  Kalpha_b_Kv43 = 12.8433;
  V0alpha_b_Kv43 = -111.33209;
  Abeta_b_Kv43 = 0.0345;
  Kbeta_b_Kv43 = -8.90123;
  V0beta_b_Kv43 = -49.9537;
  alpha_a_Kv43 = Aalpha_a_Kv43.*Q10_Kv43./(1 + exp((-V0alpha_a_Kv43 +...
    v)./Kalpha_a_Kv43));
  beta_a_Kv43 = Abeta_a_Kv43.*Q10_Kv43.*exp(-(-V0beta_a_Kv43 + v)./Kbeta_a_Kv43);
  alpha_b_Kv43 = Aalpha_b_Kv43.*Q10_Kv43./(1 + exp((-V0alpha_b_Kv43 +...
    v)./Kalpha_b_Kv43));
  beta_b_Kv43 = Abeta_b_Kv43.*Q10_Kv43./(1 + exp((-V0beta_b_Kv43 +...
    v)./Kbeta_b_Kv43));
  a_inf_Kv43 = alpha_a_Kv43./(alpha_a_Kv43 + beta_a_Kv43);
  tau_a_Kv43 = 1.0./(alpha_a_Kv43 + beta_a_Kv43);
  b_inf_Kv43 = alpha_b_Kv43./(alpha_b_Kv43 + beta_b_Kv43);
  tau_b_Kv43 = 1.0./(alpha_b_Kv43 + beta_b_Kv43);
  values(20:Ns:end) = (-a_Kv43 + a_inf_Kv43).*allow_change./tau_a_Kv43;
  values(21:Ns:end) = (-b_Kv43 + b_inf_Kv43).*allow_change./tau_b_Kv43;
  I_Kv43 = g_Kv43.*a_Kv43.^3.*(-E_K + v).*b_Kv43;

  % Expressions for the Kir2.x component
  Q10_Kir2x = 3.^(-2.0 + 0.1.*Temp);
  Aalpha_d_Kir2x = 0.13289;
  Kalpha_d_Kir2x = -24.3902;
  V0alpha_d_Kir2x = -83.94;
  Abeta_d_Kir2x = 0.16994;
  Kbeta_d_Kir2x = 35.714;
  V0beta_d_Kir2x = -83.94;
  a_d_Kir2x = Aalpha_d_Kir2x.*Q10_Kir2x.*exp((-V0alpha_d_Kir2x +...
    v)./Kalpha_d_Kir2x);
  b_d_Kir2x = Abeta_d_Kir2x.*Q10_Kir2x.*exp((-V0beta_d_Kir2x + v)./Kbeta_d_Kir2x);
  tau_d_Kir2x = 1.0./(a_d_Kir2x + b_d_Kir2x);
  d_inf_Kir2x = a_d_Kir2x./(a_d_Kir2x + b_d_Kir2x);
  values(22:Ns:end) = (-d_Kir2x + d_inf_Kir2x).*allow_change./tau_d_Kir2x;
  I_Kir2x = g_Kir2x.*(-E_K + v).*d_Kir2x;

  % Expressions for the Kca1.1 component
  q10_Kca11 = 3;
  qt_Kca11 = q10_Kca11.^(-2.3 + 0.1.*Temp);
  Qo_Kca11 = 0.73;
  Qc_Kca11 = -0.67;
  k1_Kca11 = 1000.0;
  onoffrate_Kca11 = 1;
  Kc_Kca11 = 0.011;
  Ko_Kca11 = 0.0011;
  pf0_Kca11 = 0.00239;
  pf1_Kca11 = 0.007;
  pf2_Kca11 = 0.04;
  pf3_Kca11 = 0.295;
  pf4_Kca11 = 0.557;
  pb0_Kca11 = 3.936;
  pb1_Kca11 = 1.152;
  pb2_Kca11 = 0.659;
  pb3_Kca11 = 0.486;
  pb4_Kca11 = 0.092;
  c01_Kca11 = 4.*cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c12_Kca11 = 3.*cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c23_Kca11 = 2.*cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c34_Kca11 = cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o01_Kca11 = 4.*cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o12_Kca11 = 3.*cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o23_Kca11 = 2.*cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o34_Kca11 = cai.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c10_Kca11 = Kc_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c21_Kca11 = 2.*Kc_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c32_Kca11 = 3.*Kc_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  c43_Kca11 = 4.*Kc_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o10_Kca11 = Ko_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o21_Kca11 = 2.*Ko_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o32_Kca11 = 3.*Ko_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  o43_Kca11 = 4.*Ko_Kca11.*k1_Kca11.*onoffrate_Kca11.*qt_Kca11;
  alpha_Kca11 = exp(FARADAY.*Qo_Kca11.*v./(R.*(273.15 + Temp)));
  beta_Kca11 = exp(FARADAY.*Qc_Kca11.*v./(R.*(273.15 + Temp)));
  f0_Kca11 = pf0_Kca11.*alpha_Kca11.*qt_Kca11;
  f1_Kca11 = pf1_Kca11.*alpha_Kca11.*qt_Kca11;
  f2_Kca11 = pf2_Kca11.*alpha_Kca11.*qt_Kca11;
  f3_Kca11 = pf3_Kca11.*alpha_Kca11.*qt_Kca11;
  f4_Kca11 = pf4_Kca11.*alpha_Kca11.*qt_Kca11;
  b0_Kca11 = pb0_Kca11.*beta_Kca11.*qt_Kca11;
  b1_Kca11 = pb1_Kca11.*beta_Kca11.*qt_Kca11;
  b2_Kca11 = pb2_Kca11.*beta_Kca11.*qt_Kca11;
  b3_Kca11 = pb3_Kca11.*beta_Kca11.*qt_Kca11;
  b4_Kca11 = pb4_Kca11.*beta_Kca11.*qt_Kca11;
  C4_Kca11 = 1 - C0_Kca11 - C1_Kca11 - C2_Kca11 - C3_Kca11 - O0_Kca11 -...
    O1_Kca11 - O2_Kca11 - O3_Kca11 - O4_Kca11;
  values(23:Ns:end) = ((-c01_Kca11 - f0_Kca11).*C0_Kca11 + C1_Kca11.*c10_Kca11 +...
    O0_Kca11.*b0_Kca11).*allow_change;
  values(24:Ns:end) = ((-c10_Kca11 - c12_Kca11 - f1_Kca11).*C1_Kca11 +...
    C0_Kca11.*c01_Kca11 + C2_Kca11.*c21_Kca11 + O1_Kca11.*b1_Kca11).*allow_change;
  values(25:Ns:end) = ((-c21_Kca11 - c23_Kca11 - f2_Kca11).*C2_Kca11 +...
    C1_Kca11.*c12_Kca11 + C3_Kca11.*c32_Kca11 + O2_Kca11.*b2_Kca11).*allow_change;
  values(26:Ns:end) = ((-c32_Kca11 - c34_Kca11 - f3_Kca11).*C3_Kca11 +...
    C2_Kca11.*c23_Kca11 + C4_Kca11.*c43_Kca11 + O3_Kca11.*b3_Kca11).*allow_change;
  values(27:Ns:end) = ((-b0_Kca11 - o01_Kca11).*O0_Kca11 + C0_Kca11.*f0_Kca11 +...
    O1_Kca11.*o10_Kca11).*allow_change;
  values(28:Ns:end) = ((-b1_Kca11 - o10_Kca11 - o12_Kca11).*O1_Kca11 +...
    C1_Kca11.*f1_Kca11 + O0_Kca11.*o01_Kca11 + O2_Kca11.*o21_Kca11).*allow_change;
  values(29:Ns:end) = ((-b2_Kca11 - o21_Kca11 - o23_Kca11).*O2_Kca11 +...
    C2_Kca11.*f2_Kca11 + O1_Kca11.*o12_Kca11 + O3_Kca11.*o32_Kca11).*allow_change;
  values(30:Ns:end) = ((-b3_Kca11 - o32_Kca11 - o34_Kca11).*O3_Kca11 +...
    C3_Kca11.*f3_Kca11 + O2_Kca11.*o23_Kca11 + O4_Kca11.*o43_Kca11).*allow_change;
  values(31:Ns:end) = ((-b4_Kca11 - o43_Kca11).*O4_Kca11 + C4_Kca11.*f4_Kca11 +...
    O3_Kca11.*o34_Kca11).*allow_change;
  I_Kca11 = g_Kca11.*(-E_K + v).*(O0_Kca11 + O1_Kca11 + O2_Kca11 + O3_Kca11 +...
    O4_Kca11);

  % Expressions for the Kca2.2 component
  Q10_Kca22 = 3.0;
  tcorr_Kca22 = Q10_Kca22.^(-2.3 + 0.1.*Temp);
  invc1_Kca22 = 0.08;
  invc2_Kca22 = 0.08;
  invc3_Kca22 = 0.2;
  invo1_Kca22 = 1;
  invo2_Kca22 = 0.1;
  diro1_Kca22 = 0.16;
  diro2_Kca22 = 1.2;
  dirc2_Kca22 = 200.0;
  dirc3_Kca22 = 160.0;
  dirc4_Kca22 = 80.0;
  invc1_t_Kca22 = invc1_Kca22.*tcorr_Kca22;
  invc2_t_Kca22 = invc2_Kca22.*tcorr_Kca22;
  invc3_t_Kca22 = invc3_Kca22.*tcorr_Kca22;
  invo1_t_Kca22 = invo1_Kca22.*tcorr_Kca22;
  invo2_t_Kca22 = invo2_Kca22.*tcorr_Kca22;
  diro1_t_Kca22 = diro1_Kca22.*tcorr_Kca22;
  diro2_t_Kca22 = diro2_Kca22.*tcorr_Kca22;
  dirc2_t_Kca22 = dirc2_Kca22.*tcorr_Kca22;
  dirc3_t_Kca22 = dirc3_Kca22.*tcorr_Kca22;
  dirc4_t_Kca22 = dirc4_Kca22.*tcorr_Kca22;
  dirc2_t_ca_Kca22 = cai.*dirc2_t_Kca22;
  dirc3_t_ca_Kca22 = cai.*dirc3_t_Kca22;
  dirc4_t_ca_Kca22 = cai.*dirc4_t_Kca22;
  c4_Kca22 = 1 - c1_Kca22 - c2_Kca22 - c3_Kca22 - o1_Kca22 - o2_Kca22;
  values(32:Ns:end) = (c2_Kca22.*invc1_t_Kca22 -...
    c1_Kca22.*dirc2_t_ca_Kca22).*allow_change;
  values(33:Ns:end) = ((-dirc3_t_ca_Kca22 - invc1_t_Kca22).*c2_Kca22 +...
    c1_Kca22.*dirc2_t_ca_Kca22 + c3_Kca22.*invc2_t_Kca22).*allow_change;
  values(34:Ns:end) = ((-dirc4_t_ca_Kca22 - diro1_t_Kca22 - invc2_t_Kca22).*c3_Kca22 +...
    c2_Kca22.*dirc3_t_ca_Kca22 + c4_Kca22.*invc3_t_Kca22 +...
    invo1_t_Kca22.*o1_Kca22).*allow_change;
  values(35:Ns:end) = (c3_Kca22.*diro1_t_Kca22 - invo1_t_Kca22.*o1_Kca22).*allow_change;
  values(36:Ns:end) = (c4_Kca22.*diro2_t_Kca22 - invo2_t_Kca22.*o2_Kca22).*allow_change;
  I_Kca22 = g_Kca22.*(-E_K + v).*(o1_Kca22 + o2_Kca22);

  % Expressions for the Kca3.1 component
  Yvdep_Kca31 = 13.3643751073.*exp(0.037037037037.*v);
  Yconcdep_Kca31 = ((cai < 0.01).*((7.5 - 500.*cai)./(-1 +...
    102586.491213.*exp(-769.230769231.*cai))) + ~(cai <...
    0.01).*(0.0545700593113));
  Ybeta_Kca31 = 0.05;
  Yalpha_Kca31 = Yconcdep_Kca31.*Yvdep_Kca31;
  values(37:Ns:end) = ((1 - Y_Kca31).*Yalpha_Kca31 - Ybeta_Kca31.*Y_Kca31).*allow_change;
  I_Kca31 = g_Kca31.*(-E_K + v).*Y_Kca31;

  % Expressions for the Cav2.1 component
  q10_Cav21 = 3;
  qt_Cav21 = q10_Cav21.^(-2.3 + 0.1.*Temp);
  vhalfm_Cav21 = -29.458;
  cvm_Cav21 = 8.429;
  vshift_Cav21 = 0;
  minf_Cav21 = 1.0./(1.0 + exp((vhalfm_Cav21 + vshift_Cav21 - v)./cvm_Cav21));
  taum_Cav21 = 1.0.*((v >= -40).*(0.2702 + 1.1622.*exp(0.00609050490286.*(26.798 +...
    v).*(-26.798 - v))) + ~(v >=...
    -40).*(0.6923.*exp(0.000917960072409.*v)))./qt_Cav21;
  values(38:Ns:end) = (-m_Cav21 + minf_Cav21).*allow_change./taum_Cav21;
  z_Ca = 2;
  zeta_Ca = FARADAY.*z_Ca.*v./(R.*(273.19 + Temp));
  ghk = ((abs(1 - exp(-zeta_Ca)) < 1e-06).*(1e-06.*F.*z_Ca.*(1 + zeta_Ca./2).*(cai...
    - cao.*exp(-zeta_Ca))) + ~(abs(1 - exp(-zeta_Ca)) <...
    1e-06).*(1e-06.*F.*z_Ca.*(cai - cao.*exp(-zeta_Ca)).*zeta_Ca./(1 -...
    exp(-zeta_Ca))));
  I_Cav21 = 1000.0.*g_Cav21.*m_Cav21.^3.*ghk;

  % Expressions for the Cav3.1 component
  q10_Cav31 = 3;
  qt_Cav31 = q10_Cav31.^(-3.7 + 0.1.*Temp);
  v0_m_inf_Cav31 = -52.0;
  v0_h_inf_Cav31 = -72.0;
  k_m_inf_Cav31 = -5.0;
  k_h_inf_Cav31 = 7.0;
  C_tau_m_Cav31 = 1.0;
  A_tau_m_Cav31 = 1.0;
  v0_tau_m1_Cav31 = -40.0;
  v0_tau_m2_Cav31 = -102.0;
  k_tau_m1_Cav31 = 9.0;
  k_tau_m2_Cav31 = -18.0;
  C_tau_h_Cav31 = 15.0;
  A_tau_h_Cav31 = 1.0;
  v0_tau_h1_Cav31 = -32.0;
  k_tau_h1_Cav31 = 7.0;
  minf_Cav31 = 1.0./(1 + exp((-v0_m_inf_Cav31 + v)./k_m_inf_Cav31));
  hinf_Cav31 = 1.0./(1 + exp((-v0_h_inf_Cav31 + v)./k_h_inf_Cav31));
  taum_Cav31 = ((v <= -90).*(1) + ~(v <= -90).*((C_tau_m_Cav31 +...
    A_tau_m_Cav31./(exp((-v0_tau_m1_Cav31 + v)./k_tau_m1_Cav31) +...
    exp((-v0_tau_m2_Cav31 + v)./k_tau_m2_Cav31)))./qt_Cav31));
  tauh_Cav31 = (C_tau_h_Cav31 + A_tau_h_Cav31.*exp(-(-v0_tau_h1_Cav31 +...
    v)./k_tau_h1_Cav31))./qt_Cav31;
  values(39:Ns:end) = (-m_Cav31 + minf_Cav31).*allow_change./taum_Cav31;
  values(40:Ns:end) = (-h_Cav31 + hinf_Cav31).*allow_change./tauh_Cav31;
  I_Cav31 = 1000.0.*g_Cav31.*m_Cav31.^2.*ghk.*h_Cav31;

  % Expressions for the Cav3.2 component
  phi_m_Cav32 = 6.89864830731;
  phi_h_Cav32 = 3.73719281885;
  m_inf_Cav32 = 1.0./(1 + 0.000607957605998.*exp(-0.135135135135.*v));
  h_inf_Cav32 = 1.0./(1 + 148461.062051.*exp(0.139275766017.*v));
  tau_m_Cav32 = (1.9 + 1.0./(22.4040937191.*exp(0.0840336134454.*v) +...
    0.00189854653589.*exp(-v./21)))./phi_m_Cav32;
  tau_h_Cav32 = 13.7 + (1942.0 +...
    55178666.3027.*exp(0.108695652174.*v))./(phi_h_Cav32.*(1 +...
    30321871936.2.*exp(0.27027027027.*v)));
  values(41:Ns:end) = (-m_Cav32 + m_inf_Cav32).*allow_change./tau_m_Cav32;
  values(42:Ns:end) = (-h_Cav32 + h_inf_Cav32).*allow_change./tau_h_Cav32;
  E_Ca = R.*(273.15 + Temp).*log(cao./cai)./(2.*FARADAY);
  I_Cav32 = g_Cav32.*m_Cav32.^2.*(-E_Ca + v).*h_Cav32;

  % Expressions for the Cav3.3 component
  q10_Cav33 = 2.3;
  qt_Cav33 = q10_Cav33.^(-2.8 + 0.1.*Temp);
  gCav3_3bar = 1e-05;
  vhalfn_Cav33 = -41.5;
  vhalfl_Cav33 = -69.8;
  kn_Cav33 = 6.2;
  kl_Cav33 = -6.1;
  n_inf_Cav33 = 1.0./(1.0 + exp((vhalfn_Cav33 - v)./kn_Cav33));
  l_inf_Cav33 = 1.0./(1.0 + exp((vhalfl_Cav33 - v)./kl_Cav33));
  tau_n_Cav33 = ((v > -60).*(7.2 + 0.02.*exp(-0.0680272108844.*v)) + ~(v >...
    -60).*(79.5 + 2.0.*exp(-0.10752688172.*v)))./qt_Cav33;
  tau_l_Cav33 = ((v > -60).*(0.875.*exp(120./41 + v./41)) + ~(v >...
    -60).*(260.0))./qt_Cav33;
  w_Cav33 = FARADAY.*z_Ca.*v./(R.*(273.14 + Temp));
  ghk_Cav33 = -FARADAY.*z_Ca.*(cao - cai.*exp(w_Cav33)).*w_Cav33./(-1.0 +...
    exp(w_Cav33));
  values(43:Ns:end) = (-n_Cav33 + n_inf_Cav33).*allow_change./tau_n_Cav33;
  values(44:Ns:end) = (-l_Cav33 + l_inf_Cav33).*allow_change./tau_l_Cav33;
  I_Cav33 = gCav3_3bar.*g_Cav33.*n_Cav33.^3.*ghk_Cav33.*l_Cav33;

  % Expressions for the HCN1 component
  q10_HCN1 = 3.0;
  qt_HCN1 = q10_HCN1.^(-3.7 + 0.1.*Temp);
  ratetau_HCN1 = 1.0;
  ljp_HCN1 = 9.3;
  v_inf_half_noljp_HCN1 = -90.3;
  v_inf_k_HCN1 = 9.67;
  v_tau_const_HCN1 = 0.0018;
  v_tau_half1_noljp_HCN1 = -68.0;
  v_tau_half2_noljp_HCN1 = -68.0;
  v_tau_k1_HCN1 = -22.0;
  v_tau_k2_HCN1 = 7.14;
  v_inf_half_HCN1 = v_inf_half_noljp_HCN1 - ljp_HCN1;
  v_tau_half1_HCN1 = v_tau_half1_noljp_HCN1 - ljp_HCN1;
  v_tau_half2_HCN1 = v_tau_half2_noljp_HCN1 - ljp_HCN1;
  hinf_HCN1 = 1.0./(1 + exp((-v_inf_half_HCN1 + v)./v_inf_k_HCN1));
  tauh_HCN1 = ratetau_HCN1./(v_tau_const_HCN1.*(exp((-v_tau_half1_HCN1 +...
    v)./v_tau_k1_HCN1) + exp((-v_tau_half2_HCN1 + v)./v_tau_k2_HCN1)).*qt_HCN1);
  values(45:Ns:end) = (-h_HCN1 + hinf_HCN1).*allow_change./tauh_HCN1;
  I_HCN1 = g_HCN1.*(-E_h + v).*h_HCN1;

  % Expressions for the leak component
  I_leak = g_leak.*(-E_leak + v);

  % Expressions for the Stimulation component
  stim_period = 1000.0./stim_frequency;
  I_stim = ((t - floor(t./stim_period).*stim_period <= stim_duration +...
    stim_start & t - floor(t./stim_period).*stim_period >=...
    stim_start).*(-stim_amplitude) + ~(t - floor(t./stim_period).*stim_period <=...
    stim_duration + stim_start & t - floor(t./stim_period).*stim_period >=...
    stim_start).*(0));

  % Expressions for the Membrane component
  I_tot = I_Cav21 + I_Cav31 + I_Cav32 + I_Cav33 + I_HCN1 + I_Kca11 + I_Kca22 +...
    I_Kca31 + I_Kir2x + I_Kv11 + I_Kv15 + I_Kv33 + I_Kv34 + I_Kv43 + I_Na +...
    I_leak + I_stim;
  values(46:Ns:end) = -1.0.*I_tot.*allow_change./Cm + G*v;
end
