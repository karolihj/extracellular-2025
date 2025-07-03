function [values] = Fabbri2017_rhs_vectorized(t, states, parameters)
  % Compute the right hand side of the
  % Human_SAN_Fabbri_Fantini_Wilders_Severi_2017 ODE

  % Assign states
  fCa=states(1,:); R=states(2,:); O=states(3,:); I=states(4,:); RI=states(5,:);...
    fTC=states(6,:); fTMC=states(7,:); fTMM=states(8,:); fCMi=states(9,:);...
    fCMs=states(10,:); fCQ=states(11,:); y=states(12,:); m=states(13,:);...
    h=states(14,:); dL=states(15,:); fL=states(16,:); dT=states(17,:); fT=states(18,:);...
    Cai=states(19,:); Ca_sub=states(20,:); Ca_nsr=states(21,:); Ca_jsr=states(22,:);...
    r_Kur=states(23,:); s_Kur=states(24,:); q=states(25,:); r=states(26,:);...
    paS=states(27,:); paF=states(28,:); piy=states(29,:); n=states(30,:);...
    a=states(31,:); Nai_=states(32,:); V_ode=states(33,:);

  % Assign parameters
  C=parameters(1,:); F=parameters(2,:); R_Membrane=parameters(3,:);...
    T=parameters(4,:); clamp_mode=parameters(5,:); V_holding=parameters(6,:);...
    V_test=parameters(7,:); t_holding=parameters(8,:); t_test=parameters(9,:);...
    ACh=parameters(10,:); Iso_1_uM=parameters(11,:); P_CaL=parameters(12,:);...
    Km_fCa=parameters(13,:); alpha_fCa=parameters(14,:); V_dL=parameters(15,:);...
    k_dL=parameters(16,:); k_fL=parameters(17,:); shift_fL=parameters(18,:);...
    EC50_SR=parameters(19,:); HSR=parameters(20,:); MaxSR=parameters(21,:);...
    MinSR=parameters(22,:); kiCa=parameters(23,:); kim=parameters(24,:);...
    koCa=parameters(25,:); kom=parameters(26,:); ks=parameters(27,:);...
    K_up=parameters(28,:); P_up_basal=parameters(29,:); slope_up=parameters(30,:);...
    tau_dif_Ca=parameters(31,:); tau_tr=parameters(32,:); CM_tot=parameters(33,:);...
    CQ_tot=parameters(34,:); Mgi=parameters(35,:); TC_tot=parameters(36,:);...
    TMC_tot=parameters(37,:); kb_CM=parameters(38,:); kb_CQ=parameters(39,:);...
    kb_TC=parameters(40,:); kb_TMC=parameters(41,:); kb_TMM=parameters(42,:);...
    kf_CM=parameters(43,:); kf_CQ=parameters(44,:); kf_TC=parameters(45,:);...
    kf_TMC=parameters(46,:); kf_TMM=parameters(47,:); L_cell=parameters(48,:);...
    L_sub=parameters(49,:); R_cell=parameters(50,:); V_i_part=parameters(51,:);...
    V_jsr_part=parameters(52,:); V_nsr_part=parameters(53,:); Cao=parameters(54,:);...
    Ki=parameters(55,:); Ko=parameters(56,:); Nao=parameters(57,:);...
    Nai_clamp=parameters(58,:); Km_f=parameters(59,:); alpha=parameters(60,:);...
    blockade=parameters(61,:); g_f=parameters(62,:); y_shift=parameters(63,:);...
    Km_Kp=parameters(64,:); Km_Nap=parameters(65,:); i_NaK_max=parameters(66,:);...
    K1ni=parameters(67,:); K1no=parameters(68,:); K2ni=parameters(69,:);...
    K2no=parameters(70,:); K3ni=parameters(71,:); K3no=parameters(72,:);...
    K_NaCa=parameters(73,:); Kci=parameters(74,:); Kcni=parameters(75,:);...
    Kco=parameters(76,:); Qci=parameters(77,:); Qco=parameters(78,:);...
    Qn=parameters(79,:); blockade_NaCa=parameters(80,:); g_Na=parameters(81,:);...
    g_Na_L=parameters(82,:); delta_m=parameters(83,:); P_CaT=parameters(84,:);...
    offset_fT=parameters(85,:); g_Kur=parameters(86,:); g_to=parameters(87,:);...
    g_Kr=parameters(88,:); g_Ks_=parameters(89,:); ACh_on=parameters(90,:);...
    g_KACh=parameters(91,:);

  % Init return args
  values = zeros(size(states));

  % Expressions for the Voltage clamp component
  RTONF = R_Membrane.*T./F;
  V_clamp = ((t < t_holding + t_test & t > t_holding).*(V_test) + ~(t <...
    t_holding + t_test & t > t_holding).*(V_holding));
  V = ((clamp_mode >= 1.0).*(V_clamp) + ~(clamp_mode >= 1.0).*(V_ode));

  % Expressions for the FCa gate component
  fCa_infinity = Km_fCa./(Km_fCa + Ca_sub);
  tau_fCa = 0.001.*fCa_infinity./alpha_fCa;
  values(1,:) = (-fCa + fCa_infinity)./tau_fCa;

  % Expressions for the Ca SR release component
  j_SRCarel = ks.*(-Ca_sub + Ca_jsr).*O;
  kCaSR = MaxSR - (MaxSR - MinSR)./(1.0 + (EC50_SR./Ca_jsr).^HSR);
  koSRCa = koCa./kCaSR;
  kiSRCa = kiCa.*kCaSR;
  values(2,:) = 1e-3.*(kim.*RI + kom.*O - Ca_sub.^2.0.*R.*koSRCa - Ca_sub.*R.*kiSRCa);
  values(3,:) = 1e-3.*(kim.*I - kom.*O + Ca_sub.^2.0.*R.*koSRCa - Ca_sub.*O.*kiSRCa);
  values(4,:) = 1e-3.*(-kim.*I - kom.*I + Ca_sub.^2.0.*RI.*koSRCa + Ca_sub.*O.*kiSRCa);
  values(5,:) = 1e-3.*(kom.*I - kim.*RI + Ca_sub.*R.*kiSRCa - Ca_sub.^2.0.*RI.*koSRCa);

  % Expressions for the Ca intracellular fluxes component
  b_up = ((Iso_1_uM > 0).*(-0.25) + ~(Iso_1_uM > 0).*(((ACh >...
    0).*(0.7.*ACh./(9e-05 + ACh)) + ~(ACh > 0).*(0))));
  P_up = P_up_basal.*(1.0 - b_up);
  j_Ca_dif = (-Cai + Ca_sub)./tau_dif_Ca;
  j_up = P_up./(1.0 + exp((K_up - Cai)./slope_up));
  j_tr = (-Ca_jsr + Ca_nsr)./tau_tr;

  % Expressions for the Ca buffering component
  delta_fTC = -kb_TC.*fTC + kf_TC.*(1.0 - fTC).*Cai;
  delta_fTMC = -kb_TMC.*fTMC + kf_TMC.*(1.0 - fTMC - fTMM).*Cai;
  delta_fTMM = -kb_TMM.*fTMM + Mgi.*kf_TMM.*(1.0 - fTMC - fTMM);
  delta_fCMi = -kb_CM.*fCMi + kf_CM.*(1.0 - fCMi).*Cai;
  delta_fCMs = -kb_CM.*fCMs + kf_CM.*(1.0 - fCMs).*Ca_sub;
  delta_fCQ = -kb_CQ.*fCQ + kf_CQ.*(1.0 - fCQ).*Ca_jsr;
  values(6,:) = 1e-3.*(delta_fTC);
  values(7,:) = 1e-3.*(delta_fTMC);
  values(8,:) = 1e-3.*(delta_fTMM);
  values(9,:) = 1e-3.*(delta_fCMi);
  values(10,:) = 1e-3.*(delta_fCMs);
  values(11,:) = 1e-3.*(delta_fCQ);

  % Expressions for the Cell parameters component
  V_cell = 1e-09.*pi.*L_cell.*R_cell.^2.0;
  V_sub = 2e-09.*pi.*L_cell.*L_sub.*(R_cell - 0.5.*L_sub);
  V_jsr = V_jsr_part.*V_cell;
  V_i = -V_sub + V_i_part.*V_cell;
  V_nsr = V_nsr_part.*V_cell;

  % Expressions for the Ionic values component
  Nai = Nai_;
  E_Na = RTONF.*log(Nao./Nai);
  E_K = RTONF.*log(Ko./Ki);

  % Expressions for the i_f component
  G_f = g_f./(Ko.*(Km_f + Ko));
  G_f_K = G_f./(1.0 + alpha);
  G_f_Na = alpha.*G_f_K;
  g_f_Na = Ko.*G_f_Na./(Km_f + Ko);
  g_f_K = Ko.*G_f_K./(Km_f + Ko);
  i_fNa = (1.0 - blockade).*(-E_Na + V).*g_f_Na.*y;
  i_fK = (1.0 - blockade).*(-E_K + V).*g_f_K.*y;
  i_f = i_fK + i_fNa;

  % Expressions for the y gate component
  ACh_shift = ((ACh > 0).*(-1.0 - 9.898.*ACh.^0.618./(0.00122423 +...
    1.0.*ACh.^0.618)) + ~(ACh > 0).*(0));
  Iso_shift = ((Iso_1_uM > 0).*(7.5) + ~(Iso_1_uM > 0).*(0));
  tau_y = -0.054 + 1.0./((8.73 + 0.1.*V - 0.1.*ACh_shift - 0.1.*Iso_shift)./(1.0 -...
    2.61347497529e-08.*exp(0.2.*ACh_shift + 0.2.*Iso_shift - 0.2.*V)) + (53.568 +...
    0.36.*V - 0.36.*ACh_shift - 0.36.*Iso_shift)./(-1.0 +...
    18412.7750706.*exp(0.066.*V - 0.066.*ACh_shift - 0.066.*Iso_shift)));
  y_infinity = ((V < -80.0 + y_shift + ACh_shift + Iso_shift).*(0.01329 +...
    0.99921./(1.0 + 144573.626454.*exp(0.122321166455.*V -...
    0.122321166455.*y_shift - 0.122321166455.*ACh_shift -...
    0.122321166455.*Iso_shift))) + ~(V < -80.0 + y_shift + ACh_shift +...
    Iso_shift).*(0.0002501.*exp(0.0777544514423.*y_shift +...
    0.0777544514423.*ACh_shift + 0.0777544514423.*Iso_shift -...
    0.0777544514423.*V)));
  values(12,:) = 1e-3.*((-y + y_infinity)./tau_y);

  % Expressions for the i_NaK component
  Iso_increase = ((Iso_1_uM > 0).*(1.2) + ~(Iso_1_uM > 0).*(1.0));
  i_NaK = i_NaK_max.*(1.0 + (Km_Kp./Ko).^1.2).^(-1.0).*(1.0 +...
    (Km_Nap./Nai).^1.3).^(-1.0).*(1.0 + 0.00408677143846.*exp(0.05.*E_Na -...
    0.05.*V)).^(-1.0).*Iso_increase;

  % Expressions for the i_NaCa component
  k43 = Nai./(K3ni + Nai);
  k41 = exp(-0.5.*Qn.*V./RTONF);
  di = 1.0 + (1.0 + (1.0 + Nai./K3ni).*Nai./K2ni).*Nai./K1ni + (1.0 + Nai./Kcni +...
    exp(-Qci.*V./RTONF)).*Ca_sub./Kci;
  k34 = Nao./(K3no + Nao);
  k32 = exp(0.5.*Qn.*V./RTONF);
  do_ = 1.0 + Cao.*(1.0 + exp(Qco.*V./RTONF))./Kco + Nao.*(1.0 + Nao.*(1.0 +...
    Nao./K3no)./K2no)./K1no;
  k12 = Ca_sub.*exp(-Qci.*V./RTONF)./(Kci.*di);
  k14 = Nai.^2.*(1.0 + Nai./K3ni).*exp(0.5.*Qn.*V./RTONF)./(K1ni.*K2ni.*di);
  k21 = Cao.*exp(Qco.*V./RTONF)./(Kco.*do_);
  k23 = Nao.^2.*(1.0 + Nao./K3no).*exp(-0.5.*Qn.*V./RTONF)./(K1no.*K2no.*do_);
  x1 = (k21 + k23).*k34.*k41 + (k41 + k43).*k21.*k32;
  x2 = (k12 + k14).*k32.*k43 + (k32 + k34).*k12.*k41;
  x3 = (k21 + k23).*k14.*k43 + (k41 + k43).*k12.*k23;
  x4 = (k12 + k14).*k23.*k34 + (k32 + k34).*k14.*k21;
  i_NaCa = K_NaCa.*(1.0 - blockade_NaCa).*(k21.*x2 - k12.*x1)./(x1 + x2 + x3 + x4);

  % Expressions for the i_Na component
  E_mh = RTONF.*log((Nao + 0.12.*Ko)./(0.12.*Ki + Nai));
  i_Na_ = g_Na.*m.^3.0.*(-E_mh + V).*h;
  i_Na_L = g_Na_L.*m.^3.0.*(-E_mh + V);
  i_Na = i_Na_ + i_Na_L;

  % Expressions for the m gate component
  m_infinity = 1.0./(1.0 + 0.00634650333102.*exp(-0.120328255481.*V));
  E0_m = 41.0 + V;
  alpha_m = ((abs(E0_m) < delta_m).*(2000.0) + ~(abs(E0_m) <...
    delta_m).*(200.0.*E0_m./(1.0 - exp(-0.1.*E0_m))));
  beta_m = 198.580949027.*exp(-0.056.*V);
  tau_m = 1.0./(alpha_m + beta_m);
  values(13,:) = 1e-3.*((-m + m_infinity)./tau_m);

  % Expressions for the h gate component
  h_infinity = 1.0./(1.0 + 6346493.36441.*exp(0.224391338494.*V));
  alpha_h = 0.00169636470493.*exp(-0.125.*V);
  beta_h = 2000.0./(1.0 + 0.176986998447.*exp(-0.1.*V));
  tau_h = 1.0./(alpha_h + beta_h);
  values(14,:) = 1e-3.*((-h + h_infinity)./tau_h);

  % Expressions for the i_CaL component
  Iso_increase = ((Iso_1_uM > 0).*(1.23) + ~(Iso_1_uM > 0).*(1.0));
  i_siCa = 2.0.*P_CaL.*(-Cao.*exp(-2.0.*V./RTONF) + Ca_sub).*V.*dL.*fCa.*fL./((1.0 -...
    exp(-2.0.*V./RTONF)).*RTONF);
  i_siK = 0.000365.*P_CaL.*(Ki - Ko.*exp(-1.0.*V./RTONF)).*V.*dL.*fCa.*fL./((1.0 -...
    exp(-1.0.*V./RTONF)).*RTONF);
  i_siNa = 1.85e-05.*P_CaL.*(-Nao.*exp(-1.0.*V./RTONF) + Nai).*V.*dL.*fCa.*fL./((1.0 -...
    exp(-1.0.*V./RTONF)).*RTONF);
  ACh_block = 0.31.*ACh./(9e-05 + ACh);
  i_CaL = 1.0.*(1.0 - ACh_block).*(i_siCa + i_siK + i_siNa).*Iso_increase;

  % Expressions for the DL gate component
  Iso_shift_dL = ((Iso_1_uM > 0).*(-8.0) + ~(Iso_1_uM > 0).*(0));
  Iso_slope_dL = ((Iso_1_uM > 0).*(-27.0) + ~(Iso_1_uM > 0).*(0));
  dL_infinity = 1.0./(1.0 + exp((V_dL - V + Iso_shift_dL)./(k_dL.*(1.0 +...
    0.01.*Iso_slope_dL))));
  adVm = ((V == -41.8).*(-41.80001) + ~(V == -41.8).*(((V == 0).*(0) + ~(V ==...
    0).*(((V == -6.8).*(-6.80001) + ~(V == -6.8).*(V))))));
  bdVm = ((V == -1.8).*(-1.80001) + ~(V == -1.8).*(V));
  alpha_dL = (-1.186702 - 0.02839.*adVm)./(-1.0 +...
    5.47767501694e-08.*exp(-0.4.*adVm)) - (0.57732 + 0.0849.*adVm)./(-1.0 +...
    0.242521074636.*exp(-0.208333333333.*adVm));
  beta_dL = (0.020574 + 0.01143.*bdVm)./(-1.0 + 2.05443321064.*exp(0.4.*bdVm));
  tau_dL = 0.001./(alpha_dL + beta_dL);
  values(15,:) = 1e-3.*((-dL + dL_infinity)./tau_dL);

  % Expressions for the FL gate component
  fL_infinity = 1.0./(1.0 + exp((37.4 + shift_fL + V)./(5.3 + k_fL)));
  tau_fL = 0.0443 + 0.23.*exp(-(3.6 + 0.1.*V).^2.0);
  values(16,:) = 1e-3.*((-fL + fL_infinity)./tau_fL);

  % Expressions for the i_CaT component
  i_CaT = 2.0.*P_CaT.*(-Cao.*exp(-2.0.*V./RTONF) + Ca_sub).*V.*dT.*fT./((1.0 -...
    exp(-2.0.*V./RTONF)).*RTONF);

  % Expressions for the DT gate component
  dT_infinity = 1.0./(1.0 + 0.000945651581689.*exp(-0.181818181818.*V));
  tau_dT = 0.001./(3.82842850586.*exp(0.0333333333333.*V) +...
    0.297935301196.*exp(-0.0333333333333.*V));
  values(17,:) = 1e-3.*((-dT + dT_infinity)./tau_dT);

  % Expressions for the FT gate component
  fT_infinity = 1.0./(1.0 + 5113365.83285.*exp(0.263157894737.*V));
  tau_fT = offset_fT + 1.0./(2186.53556326.*exp(0.0650195058518.*V) +...
    6.77507578516.*exp(-0.0120048019208.*V));
  values(18,:) = 1e-3.*((-fT + fT_infinity)./tau_fT);

  % Expressions for the Ca dynamics component
  values(19,:) = 1e-3.*((1.0.*V_sub.*j_Ca_dif - 1.0.*V_nsr.*j_up)./V_i - CM_tot.*delta_fCMi...
    - TC_tot.*delta_fTC - TMC_tot.*delta_fTMC);
  values(20,:) = 1e-3.*(-j_Ca_dif - CM_tot.*delta_fCMs + V_jsr.*j_SRCarel./V_sub -...
    0.5.*(-2.0.*i_NaCa + i_CaT + i_siCa)./(F.*V_sub));
  values(21,:) = 1e-3.*(-V_jsr.*j_tr./V_nsr + j_up);
  values(22,:) = 1e-3.*(-j_SRCarel - CQ_tot.*delta_fCQ + j_tr);

  % Expressions for the i_Kur component
  i_Kur = g_Kur.*(-E_K + V).*r_Kur.*s_Kur;

  % Expressions for the RKur gate component
  r_Kur_infinity = 1.0./(1.0 + 0.497741497225.*exp(-0.116279069767.*V));
  tau_r_Kur = 0.0005 + 0.009./(1.0 + 1.51689679639.*exp(0.0833333333333.*V));
  values(23,:) = 1e-3.*((-r_Kur + r_Kur_infinity)./tau_r_Kur);

  % Expressions for the SKur gate component
  s_Kur_infinity = 1.0./(1.0 + 2.11700001661.*exp(0.1.*V));
  tau_s_Kur = 3.05 + 0.59./(1.0 + 403.428793493.*exp(0.1.*V));
  values(24,:) = 1e-3.*((-s_Kur + s_Kur_infinity)./tau_s_Kur);

  % Expressions for the i_to component
  i_to = g_to.*(-E_K + V).*q.*r;

  % Expressions for the q gate component
  q_infinity = 1.0./(1.0 + 43.3467083863.*exp(0.0769230769231.*V));
  tau_q = 0.00606 + 0.039102./(0.0168716780457.*exp(-0.08.*V) +...
    6.42137321286.*exp(0.1.*V));
  values(25,:) = 1e-3.*((-q + q_infinity)./tau_q);

  % Expressions for the r gate component
  r_infinity = 1.0./(1.0 + 3.62069742698.*exp(-0.0666666666667.*V));
  tau_r = 0.00275352 + 0.01440516./(16.3010892258.*exp(0.09.*V) +...
    0.0211152735604.*exp(-0.12.*V));
  values(26,:) = 1e-3.*((-r + r_infinity)./tau_r);

  % Expressions for the i_Kr component
  i_Kr = g_Kr.*(-E_K + V).*(0.1.*paS + 0.9.*paF).*piy;

  % Expressions for the Pa gate component
  pa_infinity = 1.0./(1.0 + 0.270564851304.*exp(-0.130536373961.*V));
  tau_paS = 0.84655354./(4.2.*exp(0.0588235294118.*V) +...
    0.15.*exp(-0.0462962962963.*V));
  tau_paF = 1.0./(30.0.*exp(0.1.*V) + exp(-0.0833333333333.*V));
  values(27,:) = 1e-3.*((-paS + pa_infinity)./tau_paS);
  values(28,:) = 1e-3.*((-paF + pa_infinity)./tau_paF);

  % Expressions for the Pi gate component
  tau_pi = 1.0./(100.0.*exp(-0.0182999359502.*V) + 656.0.*exp(0.00942000998521.*V));
  pi_infinity = 1.0./(1.0 + 5.32554268928.*exp(0.0584795321637.*V));
  values(29,:) = 1e-3.*((-piy + pi_infinity)./tau_pi);

  % Expressions for the i_Ks component
  g_Ks = ((Iso_1_uM > 0).*(1.2.*g_Ks_) + ~(Iso_1_uM > 0).*(g_Ks_));
  E_Ks = RTONF.*log((Ko + 0.12.*Nao)./(Ki + 0.12.*Nai));
  i_Ks = n.^2.0.*(-E_Ks + V).*g_Ks;

  % Expressions for the n gate component
  Iso_shift = ((Iso_1_uM > 0).*(-14.0) + ~(Iso_1_uM > 0).*(0));
  n_infinity = 1.0./sqrt(1.0 + 0.942127514153.*exp(0.0933959708978.*Iso_shift -...
    0.0933959708978.*V));
  alpha_n = 28.0./(1.0 + 617437.626912.*exp(0.333333333333.*Iso_shift -...
    0.333333333333.*V));
  beta_n = 1.22140275816.*exp(0.04.*Iso_shift - 0.04.*V);
  tau_n = 1.0./(alpha_n + beta_n);
  values(30,:) = 1e-3.*((-n + n_infinity)./tau_n);

  % Expressions for the i_KACh component
  i_KACh = ((ACh > 0).*(ACh_on.*g_KACh.*(1.0 + 2.71828182846.*exp(0.05.*V)).*(-E_K...
    + V).*a) + ~(ACh > 0).*(0));

  % Expressions for the a gate component
  alpha_a = 0.025641 + 3.573159./(1.0 + 1.2155e-06.*ACh.^(-1.6951));
  beta_a = 17.0233357337.*exp(0.0133.*V);
  a_infinity = alpha_a./(alpha_a + beta_a);
  tau_a = 1.0./(alpha_a + beta_a);
  values(31,:) = 1e-3.*((-a + a_infinity)./tau_a);

  % Expressions for the Nai_concentration component
  values(32,:) = 1e-3.*((-1.0 + 1.0.*Nai_clamp).*(3.0.*i_NaCa + 3.0.*i_NaK + i_Na + i_fNa...
    + i_siNa)./(F.*(1.0.*V_i + 1.0.*V_sub)));

  % Expressions for the Membrane component
  i_tot = i_CaL + i_CaT + i_KACh + i_Kr + i_Ks + i_Kur + i_Na + i_NaCa +...
    i_NaK + i_f + i_to;
  values(33,:) = 1e-3.*(-i_tot./C);
end