function [values] = base_model_rhs_vectorized(t, states, parameters)
  % Compute the right hand side of the base_model ODE

  % Assign states
  m=states(1,:); j=states(2,:); mL=states(3,:); hL=states(4,:); Xr1=states(5,:);...
    Xr2=states(6,:); x_Ks=states(7,:); q=states(8,:); r=states(9,:); d=states(10,:);...
    f=states(11,:); f_Ca_B=states(12,:); xf=states(13,:); r_RyR=states(14,:);...
    bt=states(15,:); R_mech=states(16,:); D_mech=states(17,:); X_mech=states(18,:);...
    P_mech=states(19,:); l=states(20,:); lv=states(21,:); cn=states(22,:);...
    cc=states(23,:); cd=states(24,:); csl=states(25,:); cs=states(26,:);...
    bc=states(27,:); bd=states(28,:); bs=states(29,:); bsl=states(30,:);...
    V_m=states(31,:); Na_i=states(32,:);

  % Assign parameters
  g_Na=parameters(1,:);  g_NaL=parameters(2,:);thL=parameters(3,:); ...
      KmKo=parameters(4,:); KmNaip=parameters(5,:); g_NaK=parameters(6,:); ...
      g_Kr=parameters(7,:); g_Ks=parameters(8,:); g_to=parameters(9,:); ...
      g_K1=parameters(10,:); g_bCl=parameters(11,:); g_CaL=parameters(12,:); ...
      Kdact=parameters(13,:); KmCai=parameters(14,:); KmCao=parameters(15,:); ...
      KmNai=parameters(16,:); KmNao=parameters(17,:); g_NaCa=parameters(18,:);...
    ksat=parameters(19,:); nu=parameters(20,:); KmPCa=parameters(21,:);...
    g_pCa=parameters(22,:); g_bCa=parameters(23,:);...
    E_f=parameters(24,:); g_f=parameters(25,:); Nao=parameters(26,:);...
    K_i=parameters(27,:); Ko=parameters(28,:); ce=parameters(29,:);...
    K_RyR=parameters(30,:); alpha_RyR=parameters(31,:); beta_RyR=parameters(32,:);...
    eta_RyR=parameters(33,:); gamma_RyR=parameters(34,:);...
    lambda_RyR=parameters(35,:); Vc=parameters(36,:); Vd=parameters(37,:);...
    Vn=parameters(38,:); Vs=parameters(39,:); Vsl=parameters(40,:);...
    J_SERCA_bar=parameters(41,:); K_c=parameters(42,:); K_n=parameters(43,:);...
    B_tot_c=parameters(44,:); B_tot_d=parameters(45,:); B_tot_s=parameters(46,:);...
    B_tot_sl=parameters(47,:); k_off_c=parameters(48,:); k_off_d=parameters(49,:);...
    k_off_s=parameters(50,:); k_off_sl=parameters(51,:); k_on_c=parameters(52,:);...
    k_on_d=parameters(53,:); k_on_s=parameters(54,:); k_on_sl=parameters(55,:);...
    alpha_d_c=parameters(56,:); alpha_n_s=parameters(57,:);...
    alpha_sl_c=parameters(58,:); ...
    Cli=parameters(59,:); Clo=parameters(60,:); Cm=parameters(61,:);...
    Frdy=parameters(62,:); R=parameters(63,:); Temp=parameters(64,:);...
    chi=parameters(65,:); ...
    stim_amplitude=parameters(66,:); stim_duration=parameters(67,:);...
    stim_period=parameters(68,:); stim_start=parameters(69,:); Fa=parameters(70,:);...
    H_mech=parameters(71,:); K_mech=parameters(72,:); Trop_frac=parameters(73,:);...
    a_mech=parameters(74,:); b_mech=parameters(75,:); kDX0=parameters(76,:);...
    kPD0=parameters(77,:); kPX0=parameters(78,:); kRD0=parameters(79,:);...
    kXD0=parameters(80,:); kXP0=parameters(81,:); k_off_t=parameters(82,:);...
    k_on_t=parameters(83,:); kv_mech=parameters(84,:); l0=parameters(85,:);...
    m_mech=parameters(86,:);

  % Init return args
  values = zeros(size(states));

  % Expressions for the Reversal potentials component
  FoRT = Frdy./(R.*Temp);
  ena = log(Nao./Na_i)./FoRT;
  ek = log(Ko./K_i)./FoRT;
  eca_sl = log(ce./csl)./(2.*FoRT);
  ecl = log(Cli./Clo)./FoRT;

  % Expressions for the I_Na component
  mss = (1 + exp(-19./3 - V_m./9)).^(-2);
  taum = 0.06.*exp(-(-5./51 + V_m./51).^2) + 0.13.*exp(-(23./8 + V_m./16).^2);
  aj = ((V_m >= -40).*(0) + ~(V_m >= -40).*((38 + V_m).*(-25000.0.*exp(0.2.*V_m) -...
    7e-06.*exp(-0.04.*V_m))./(1 + 19623624323.7.*exp(0.3.*V_m))));
  bj = ((V_m >= -40).*(0.6.*exp(0.09.*V_m)./(1 + exp(-40 - V_m))) + ~(V_m >=...
    -40).*(0.02.*exp(-0.01.*V_m)./(1 + 0.00369786371648.*exp(-0.14.*V_m))));
  tauj = 1.0./(aj + bj);
  jss = (1 + exp(72./7 + V_m./7)).^(-2);
  values(1,:) = (-m + mss)./taum;
  values(2,:) = (-j + jss)./tauj;
  I_Na = g_Na.*m.^3.*(-ena + V_m).*j;

  % Expressions for the I_NaL component
  mLss = 1.0./(1 + exp(-43./5 - V_m./5));
  tm = 1.0./(8.6.*exp(-77./6 - V_m./6) + 6.8.*exp(12./35 + V_m./35));
  tmL = tm;
  values(3,:) = (-mL + mLss)./tmL;
  hLss = 1.0./(1 + 124658.506952.*exp(0.133333333333.*V_m));
  values(4,:) = (-hL + hLss)./thL;
  I_NaL = (-ena + V_m).*g_NaL.*hL.*mL;

  % Expressions for the I_NaK component
  sigma = -1./7 + exp(Nao./67)./7;
  fNaK = 1.0./(1 + 0.12.*exp(-0.1.*FoRT.*V_m) + 0.037.*exp(-FoRT.*V_m).*sigma);
  I_NaK = Ko.*g_NaK.*fNaK./((1 + KmNaip.^4./Na_i.^4).*(KmKo + Ko));

  % Expressions for the I_Kr component
  Xr1_inf = 1.0./(1.0 + 0.0146327985189.*exp(-0.204081632653.*V_m));
  alpha_Xr1 = 450.0./(1.0 + 0.0111089965382.*exp(-0.1.*V_m));
  beta_Xr1 = 6.0./(1.0 + 13.5813245226.*exp(0.0869565217391.*V_m));
  tau_Xr1 = 1.0.*alpha_Xr1.*beta_Xr1;
  values(5,:) = (-Xr1 + Xr1_inf)./tau_Xr1;
  Xr2_infinity = 1.0./(1.0 + 5.8124373944.*exp(0.02.*V_m));
  alpha_Xr2 = 3.0./(1.0 + 0.0497870683679.*exp(-0.05.*V_m));
  beta_Xr2 = 1.12./(1.0 + 0.0497870683679.*exp(0.05.*V_m));
  tau_Xr2 = 1.0.*alpha_Xr2.*beta_Xr2;
  values(6,:) = (-Xr2 + Xr2_infinity)./tau_Xr2;
  I_Kr = 0.430331482912.*g_Kr.*sqrt(Ko).*(-ek + V_m).*Xr1.*Xr2;

  % Expressions for the I_Ks component
  eks = log((Ko + Nao.*0.018)./(K_i + 0.018.*Na_i))./FoRT;
  xsss = 1.0./(1 + 0.76228973079.*exp(-V_m./14));
  tauxs = 990./(1 + 0.842460441617.*exp(-V_m./14));
  values(7,:) = (-x_Ks + xsss)./tauxs;
  I_Ks = g_Ks.*x_Ks.^2.*(-eks + V_m);

  % Expressions for the i_to component
  q_inf = 1.0./(1.0 + 58.9637634804.*exp(0.0769230769231.*V_m));
  tau_q = 6 + 39.0./(0.0168716780457.*exp(-0.08.*V_m) +...
    6.46648051673.*exp(0.1.*V_m));
  values(8,:) = (-q + q_inf)./tau_q;
  r_inf = 1.0./(1.0 + 3.28489055021.*exp(-0.0533333333333.*V_m));
  tau_r = 2.75 + 14.4./(0.0207698622486.*exp(-0.12.*V_m) +...
    15.7194688773.*exp(0.09.*V_m));
  values(9,:) = (-r + r_inf)./tau_r;
  I_to = g_to.*(-ek + V_m).*q.*r;

  % Expressions for the I_K1 component
  aK1 = 1.0./(1 + 7.50455791508e-06.*exp(0.2.*V_m - 0.2.*ek));
  bK1 = (0.745912348821.*exp(0.08.*V_m - 0.08.*ek) +...
    3.32464030033e-16.*exp(0.06.*V_m - 0.06.*ek))./(1 +...
    0.0820849986239.*exp(0.5.*ek - 0.5.*V_m));
  K1ss = aK1./(aK1 + bK1);
  I_K1 = 0.430331482912.*g_K1.*sqrt(Ko).*(-ek + V_m).*K1ss;

  % Expressions for the I_bCl component
  I_bCl = g_bCl.*(-ecl + V_m);

  % Expressions for the I_Ca component
  fss = 1.0./(1 + 48.856571275.*exp(0.111111111111.*V_m)) + 0.6./(1 +...
    12.1824939607.*exp(-0.05.*V_m));
  dss = 1.0./(1 + 0.434598208507.*exp(-0.166666666667.*V_m));
  taud = (1 - 0.434598208507.*exp(-0.166666666667.*V_m)).*dss./(0.175 + 0.035.*V_m);
  tauf = 1.0./(0.02 + 0.02.*exp(-(0.493 + 0.034.*V_m).^2));
  values(10,:) = (-d + dss)./taud;
  values(11,:) = (-f + fss)./tauf;
  values(12,:) = -0.012.*f_Ca_B + 1.7.*(1 - f_Ca_B).*cd;
  ibarca_j = 4.*Frdy.*g_CaL.*(-0.34.*ce + 0.34.*cd.*exp(2.*FoRT.*V_m)).*FoRT.*V_m./(-1 +...
    exp(2.*FoRT.*V_m));
  I_CaL = (1 - f_Ca_B).*d.*f.*ibarca_j;

  % Expressions for the I_NCX component
  Ka_sl = 1.0./(1.0 + Kdact.^2./csl.^2);
  s1_sl = ce.*Na_i.^3.*exp(nu.*FoRT.*V_m);
  s2_sl = Nao.^3.*csl.*exp((-1 + nu).*FoRT.*V_m);
  s3_sl = KmCao.*Na_i.^3 + ce.*Na_i.^3 + Nao.^3.*csl + KmCai.*Nao.^3.*(1 +...
    Na_i.^3./KmNai.^3) + KmNao.^3.*(1 + csl./KmCai).*csl;
  I_NaCa = g_NaCa.*(-s2_sl + s1_sl).*Ka_sl./((1 +...
    ksat.*exp((-1 + nu).*FoRT.*V_m)).*s3_sl);

  % Expressions for the I_PCa component
  I_pCa = g_pCa.*csl.^2./(KmPCa.^2 + csl.^2);

  % Expressions for the I_CaBK component
  I_bCa = g_bCa.*(-eca_sl + V_m);

  % Expressions for the I_f component
  xf_inf = 1.0./(1 + 5956538.01318.*exp(0.2.*V_m));
  tau_xf = 1900.0./(1 + 4.48168907034.*exp(0.1.*V_m));
  values(13,:) = (-xf + xf_inf)./tau_xf;
  I_f = g_f.*(-E_f + V_m).*xf;

  % Expressions for the Ca Fluxes component
  J_CaL = -Cm.*chi.*I_CaL./(2.*Frdy);
  J_pCa = -Cm.*chi.*I_pCa./(2.*Frdy);
  J_bCa = -Cm.*chi.*I_bCa./(2.*Frdy);
  J_NaCa = Cm.*chi.*I_NaCa./Frdy;
  J_e_sl = J_NaCa + J_bCa + J_pCa;
  J_SERCA = J_SERCA_bar.*(cc.^2./K_c.^2 - cn.^2./K_n.^2)./(1 + cc.^2./K_c.^2 +...
    cn.^2./K_n.^2);
  J_n_s = alpha_n_s.*(-cs + cn);
  J_sl_c = alpha_sl_c.*(-cc + csl);
  J_d_c = alpha_d_c.*(-cc + cd);

  % Expressions for the RyRs component
  p = 1.0./(1 + K_RyR.^3./cd.^3);
  J_RyR_active = alpha_RyR.*lambda_RyR.*(-csl + cs).*p.*r_RyR;
  J_leak = alpha_RyR.*gamma_RyR.*lambda_RyR.*(-csl + cs);
  values(14,:) = eta_RyR.*(1 - r_RyR)./p -...
    J_RyR_active./(beta_RyR.*lambda_RyR);
  J_RyR = J_RyR_active + J_leak;

  % Expressions for the Mechanics component
  B_tot_t = B_tot_c.*Trop_frac;
  B_tot_c_minus_trop = B_tot_c.*(1 - Trop_frac);
  J_t_b = -k_off_t.*bt + k_on_t.*(-bt + B_tot_t).*cc;
  values(15,:) = J_t_b;
  bbt = bt./B_tot_t;
  CC = bbt.^H_mech./(K_mech.^H_mech + bbt.^H_mech);
  kRD = kRD0.*CC;
  kDX = kDX0;
  kXP = kXP0;
  kXD = kXD0;
  kPX = kPX0;
  kPD = kPD0;
  kDR = ((1.0./kRD < 100).*(1.0./kRD) + ~(1.0./kRD < 100).*(100));
  values(16,:) = D_mech.*kDR - R_mech.*kRD;
  values(17,:) = P_mech.*kPD + R_mech.*kRD + X_mech.*kXD - (kDR + kDX).*D_mech;
  values(18,:) = D_mech.*kDX + P_mech.*kPX - (kXD + kXP).*X_mech;
  values(19,:) = X_mech.*kXP - (kPD + kPX).*P_mech;
  Fp = a_mech.*((l < l0).*(1 - exp(b_mech.*(l0 - l))) + ~(l < l0).*(-1 +...
    exp(b_mech.*(-l0 + l))));
  values(20,:) = (-kv_mech.*(-l0 + l) + lv)./m_mech;
  values(21,:) = -Fp - Fa.*P_mech;

  % Expressions for the Ca Buffers component
  J_c_b = Vc.*(-k_off_c.*bc + k_on_c.*B_tot_c_minus_trop.*cc);
  J_d_b = Vd.*(-k_off_d.*bd + k_on_d.*(-bd + B_tot_d).*cd);
  J_s_b = Vs.*(-k_off_s.*bs + k_on_s.*(-bs + B_tot_s).*cs);
  J_sl_b = Vsl.*(-k_off_sl.*bsl + k_on_sl.*(-bsl +...
    B_tot_sl).*csl);

  % Expressions for the Ca Concentrations component
  values(22,:) = 1.0.*(-J_n_s + J_SERCA)./Vn;
  values(23,:) = 1.0.*(-J_SERCA - J_c_b - Vc.*J_t_b + J_d_c + J_sl_c)./Vc;
  values(24,:) = 1.0.*(-J_d_b - J_d_c + J_CaL)./Vd;
  values(25,:) = 1.0.*(-J_sl_b - J_sl_c + J_RyR + J_e_sl)./Vsl;
  values(26,:) = 1.0.*(-J_RyR - J_s_b + J_n_s)./Vs;

  % Expressions for the Ca Buffer Concentrations component
  values(27,:) = 1.0.*J_c_b./Vc;
  values(28,:) = 1.0.*J_d_b./Vd;
  values(29,:) = 1.0.*J_s_b./Vs;
  values(30,:) = 1.0.*J_sl_b./Vsl;

  % Expressions for the Membrane potential component
  % i_Stim = ((V_m < -40).*(1) + ~(V_m < -40).*(0)).*((t -...
  %   stim_period.*floor(t./stim_period) <= stim_duration + stim_start & t -...
  %   stim_period.*floor(t./stim_period) >= stim_start).*(-stim_amplitude) + ~(t -...
  %   stim_period.*floor(t./stim_period) <= stim_duration + stim_start & t -...
  %   stim_period.*floor(t./stim_period) >= stim_start).*(0));
  i_Stim = ((t -...
    stim_period.*floor(t./stim_period) <= stim_duration + stim_start & t -...
    stim_period.*floor(t./stim_period) >= stim_start).*(-stim_amplitude) + ~(t -...
    stim_period.*floor(t./stim_period) <= stim_duration + stim_start & t -...
    stim_period.*floor(t./stim_period) >= stim_start).*(0));
  I_tot = I_CaL + I_K1 + I_Kr + I_Ks + I_Na + I_NaCa + I_NaK + I_NaL + I_bCa...
    + I_bCl + I_f + I_pCa + I_to;
  values(31,:) = -I_tot - i_Stim;

  % Expressions for the Na Concentrations component
  I_Na_tot = 3.*I_NaCa + 3.*I_NaK + 0.3293.*I_f + I_Na + I_NaL;
  J_Na = -Cm.*chi.*I_Na_tot./Frdy;
  values(32,:) = J_Na;
end
