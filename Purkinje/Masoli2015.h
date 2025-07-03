#include <cmath>
#include <cstring>
#include <stdexcept>
// Gotran generated C/C++ code for the "Masoli2015" model

// Init state values
void init_state_values(double* states)
{
  states[0] = 0.745657349519; // C1_Na;
  states[1] = 0.18681749783; // C2_Na;
  states[2] = 0.0175508753271; // C3_Na;
  states[3] = 0.00073259055439; // C4_Na;
  states[4] = 1.136427441e-05; // C5_Na;
  states[5] = 4.176150917e-05; // O_Na;
  states[6] = 0.00012203268796; // B_Na;
  states[7] = 0.0054257348267; // I1_Na;
  states[8] = 0.0150081485653; // I2_Na;
  states[9] = 0.0155777872341; // I3_Na;
  states[10] = 0.00718538355943; // I4_Na;
  states[11] = 0.00123878700626; // I5_Na;
  states[12] = 0.0784327349458; // n_Kv11;
  states[13] = 0.062693414843; // m_Kv15;
  states[14] = 0.953843893364; // n_Kv15;
  states[15] = 0.982440214681; // u_Kv15;
  states[16] = 0.0206633026411; // n_Kv33;
  states[17] = 0.143344433759; // m_Kv34;
  states[18] = 0.882087088068; // h_Kv34;
  states[19] = 0.0299743553258; // a_Kv43;
  states[20] = 0.185952387828; // b_Kv43;
  states[21] = 0.206778862228; // d_Kir2x;
  states[22] = 0.983767260817; // C0_Kca11;
  states[23] = 0.0161123619009; // C1_Kca11;
  states[24] = 9.895934853e-05; // C2_Kca11;
  states[25] = 2.7012724e-07; // C3_Kca11;
  states[26] = 1.801623e-05; // O0_Kca11;
  states[27] = 2.94586477e-06; // O1_Kca11;
  states[28] = 1.8046753e-07; // O2_Kca11;
  states[29] = 4.91765e-09; // O3_Kca11;
  states[30] = 5.012e-11; // O4_Kca11;
  states[31] = 0.88748612049; // c1_Kca22;
  states[32] = 0.0999548330279; // c2_Kca22;
  states[33] = 0.00900749464037; // c3_Kca22;
  states[34] = 0.00144123766793; // o1_Kca22;
  states[35] = 0.00194798850465; // o2_Kca22;
  states[36] = 0.00342250393833; // Y_Kca31;
  states[37] = 0.0113511215154; // m_Cav21;
  states[38] = 0.176172416128; // m_Cav31;
  states[39] = 0.191716948331; // h_Cav31;
  states[40] = 0.15549099887; // m_Cav32;
  states[41] = 0.0574480599015; // h_Cav32;
  states[42] = 0.440000580693; // n_Cav33;
  states[43] = 0.164229960931; // l_Cav33;
  states[44] = 0.00144283978473; // h_HCN1;
  states[45] = -72.9572879674; // v;
}

// Default parameter values
void init_parameters_values(double* parameters)
{
  parameters[0] = 0.001; // Cm;
  parameters[1] = 96485.3365; // F;
  parameters[2] = 96.4853365; // FARADAY;
  parameters[3] = 8.31446261815; // R;
  parameters[4] = 0.0; // Tdiff;
  parameters[5] = 37; // Temp;
  parameters[6] = 4.5e-05; // cai;
  parameters[7] = 2.0; // cao;
  parameters[8] = 60; // E_Na;
  parameters[9] = 0.214; // g_Na;
  parameters[10] = -88; // E_K;
  parameters[11] = 0.002; // g_Kv11;
  parameters[12] = 0; // g_Kv15;
  parameters[13] = 0; // g_Kv33;
  parameters[14] = 0.05; // g_Kv34;
  parameters[15] = 0; // g_Kv43;
  parameters[16] = 3e-05; // g_Kir2x;
  parameters[17] = 0.01; // g_Kca11;
  parameters[18] = 0.001; // g_Kca22;
  parameters[19] = 0.01; // g_Kca31;
  parameters[20] = 0.00022; // g_Cav21;
  parameters[21] = 7e-06; // g_Cav31;
  parameters[22] = 0.0008; // g_Cav32;
  parameters[23] = 0.0001; // g_Cav33;
  parameters[24] = -34.4; // E_h;
  parameters[25] = 0.0004; // g_HCN1;
  parameters[26] = -63; // E_leak;
  parameters[27] = 0.0011; // g_leak;
}

// State index
int state_index(const char name[])
{
  // State names
  char names[][9] = {"C1_Na", "C2_Na", "C3_Na", "C4_Na", "C5_Na", "O_Na",
    "B_Na", "I1_Na", "I2_Na", "I3_Na", "I4_Na", "I5_Na", "n_Kv11", "m_Kv15",
    "n_Kv15", "u_Kv15", "n_Kv33", "m_Kv34", "h_Kv34", "a_Kv43", "b_Kv43",
    "d_Kir2x", "C0_Kca11", "C1_Kca11", "C2_Kca11", "C3_Kca11", "O0_Kca11",
    "O1_Kca11", "O2_Kca11", "O3_Kca11", "O4_Kca11", "c1_Kca22", "c2_Kca22",
    "c3_Kca22", "o1_Kca22", "o2_Kca22", "Y_Kca31", "m_Cav21", "m_Cav31",
    "h_Cav31", "m_Cav32", "h_Cav32", "n_Cav33", "l_Cav33", "h_HCN1", "v"};

  int i;
  for (i=0; i<46; i++)
  {
    if (strcmp(names[i], name)==0)
    {
      return i;
    }
  }
  return -1;
}

// Parameter index
int parameter_index(const char name[])
{
  // Parameter names
  char names[][8] = {"Cm", "F", "FARADAY", "R", "Tdiff", "Temp", "cai",
    "cao", "E_Na", "g_Na", "E_K", "g_Kv11", "g_Kv15", "g_Kv33", "g_Kv34",
    "g_Kv43", "g_Kir2x", "g_Kca11", "g_Kca22", "g_Kca31", "g_Cav21",
    "g_Cav31", "g_Cav32", "g_Cav33", "E_h", "g_HCN1", "E_leak", "g_leak"};

  int i;
  for (i=0; i<28; i++)
  {
    if (strcmp(names[i], name)==0)
    {
      return i;
    }
  }
  return -1;
}

// Compute the right hand side of the Masoli2015 ODE
void rhs(const double* states, const double t, const double* parameters,
  double* values)
{

  // Assign states
  const double C1_Na = states[0];
  const double C2_Na = states[1];
  const double C3_Na = states[2];
  const double C4_Na = states[3];
  const double C5_Na = states[4];
  const double O_Na = states[5];
  const double B_Na = states[6];
  const double I1_Na = states[7];
  const double I2_Na = states[8];
  const double I3_Na = states[9];
  const double I4_Na = states[10];
  const double I5_Na = states[11];
  const double n_Kv11 = states[12];
  const double m_Kv15 = states[13];
  const double n_Kv15 = states[14];
  const double u_Kv15 = states[15];
  const double n_Kv33 = states[16];
  const double m_Kv34 = states[17];
  const double h_Kv34 = states[18];
  const double a_Kv43 = states[19];
  const double b_Kv43 = states[20];
  const double d_Kir2x = states[21];
  const double C0_Kca11 = states[22];
  const double C1_Kca11 = states[23];
  const double C2_Kca11 = states[24];
  const double C3_Kca11 = states[25];
  const double O0_Kca11 = states[26];
  const double O1_Kca11 = states[27];
  const double O2_Kca11 = states[28];
  const double O3_Kca11 = states[29];
  const double O4_Kca11 = states[30];
  const double c1_Kca22 = states[31];
  const double c2_Kca22 = states[32];
  const double c3_Kca22 = states[33];
  const double o1_Kca22 = states[34];
  const double o2_Kca22 = states[35];
  const double Y_Kca31 = states[36];
  const double m_Cav21 = states[37];
  const double m_Cav31 = states[38];
  const double h_Cav31 = states[39];
  const double m_Cav32 = states[40];
  const double h_Cav32 = states[41];
  const double n_Cav33 = states[42];
  const double l_Cav33 = states[43];
  const double h_HCN1 = states[44];
  const double v = states[45];

  // Assign parameters
  const double Cm = parameters[0];
  const double F = parameters[1];
  const double FARADAY = parameters[2];
  const double R = parameters[3];
  const double Tdiff = parameters[4];
  const double Temp = parameters[5];
  const double cai = parameters[6];
  const double cao = parameters[7];
  const double E_Na = parameters[8];
  const double g_Na = parameters[9];
  const double E_K = parameters[10];
  const double g_Kv11 = parameters[11];
  const double g_Kv15 = parameters[12];
  const double g_Kv33 = parameters[13];
  const double g_Kv34 = parameters[14];
  const double g_Kv43 = parameters[15];
  const double g_Kir2x = parameters[16];
  const double g_Kca11 = parameters[17];
  const double g_Kca22 = parameters[18];
  const double g_Kca31 = parameters[19];
  const double g_Cav21 = parameters[20];
  const double g_Cav31 = parameters[21];
  const double g_Cav32 = parameters[22];
  const double g_Cav33 = parameters[23];
  const double E_h = parameters[24];
  const double g_HCN1 = parameters[25];
  const double E_leak = parameters[26];
  const double g_leak = parameters[27];

  // Expressions for the Nav1.6 component
  const double allow_change = (t >= Tdiff ? 1. : 0.);
  const double q10_Na = 3.;
  const double qt_Na = std::pow(q10_Na, -2.2 + 0.1*Temp);
  const double Con_Na = 0.00500000000000000;
  const double Coff_Na = 0.500000000000000;
  const double Oon_Na = 0.750000000000000;
  const double Ooff_Na = 0.00500000000000000;
  const double alpha_Na = 150.;
  const double beta_Na = 3.;
  const double gamma_Na = 150.;
  const double delta_Na = 40.;
  const double epsilon_Na = 1.75000000000000;
  const double zeta_Na = 0.0300000000000000;
  const double x1_Na = 20.;
  const double x2_Na = -20.;
  const double x3_Na = 1000000000000.00;
  const double x4_Na = -1000000000000.00;
  const double x5_Na = 1000000000000.00;
  const double x6_Na = -25.;
  const double alfac_Na = std::pow(Oon_Na/Con_Na, 0.25);
  const double btfac_Na = std::pow(Ooff_Na/Coff_Na, 0.25);
  const double f01_Na = 4.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f02_Na = 3.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f03_Na = 2.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f04_Na = alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f0O_Na = gamma_Na*std::exp(v/x3_Na)*qt_Na;
  const double fip_Na = epsilon_Na*std::exp(v/x5_Na)*qt_Na;
  const double f11_Na = 4.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f12_Na = 3.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f13_Na = 2.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f14_Na = alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f1n_Na = gamma_Na*std::exp(v/x3_Na)*qt_Na;
  const double fi1_Na = Con_Na*qt_Na;
  const double fi2_Na = Con_Na*alfac_Na*qt_Na;
  const double fi3_Na = Con_Na*(alfac_Na*alfac_Na)*qt_Na;
  const double fi4_Na = Con_Na*(alfac_Na*alfac_Na*alfac_Na)*qt_Na;
  const double fi5_Na = Con_Na*std::pow(alfac_Na, 4.)*qt_Na;
  const double fin_Na = Oon_Na*qt_Na;
  const double b01_Na = beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b02_Na = 2.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b03_Na = 3.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b04_Na = 4.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b0O_Na = delta_Na*std::exp(v/x4_Na)*qt_Na;
  const double bip_Na = zeta_Na*std::exp(v/x6_Na)*qt_Na;
  const double b11_Na = beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b12_Na = 2.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b13_Na = 3.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b14_Na = 4.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b1n_Na = delta_Na*std::exp(v/x4_Na)*qt_Na;
  const double bi1_Na = Coff_Na*qt_Na;
  const double bi2_Na = Coff_Na*btfac_Na*qt_Na;
  const double bi3_Na = Coff_Na*(btfac_Na*btfac_Na)*qt_Na;
  const double bi4_Na = Coff_Na*(btfac_Na*btfac_Na*btfac_Na)*qt_Na;
  const double bi5_Na = Coff_Na*std::pow(btfac_Na, 4.)*qt_Na;
  const double bin_Na = Ooff_Na*qt_Na;
  const double I6_Na = 1. - B_Na - C1_Na - C2_Na - C3_Na - C4_Na - C5_Na -
    I1_Na - I2_Na - I3_Na - I4_Na - I5_Na - O_Na;
  values[0] = ((-f01_Na - fi1_Na)*C1_Na + C2_Na*b01_Na +
    I1_Na*bi1_Na)*allow_change;
  values[1] = ((-b01_Na - f02_Na - fi2_Na)*C2_Na + C1_Na*f01_Na +
    C3_Na*b02_Na + I2_Na*bi2_Na)*allow_change;
  values[2] = ((-b02_Na - f03_Na - fi3_Na)*C3_Na + C2_Na*f02_Na +
    C4_Na*b03_Na + I3_Na*bi3_Na)*allow_change;
  values[3] = ((-b03_Na - f04_Na - fi4_Na)*C4_Na + C3_Na*f03_Na +
    C5_Na*b04_Na + I4_Na*bi4_Na)*allow_change;
  values[4] = ((-b04_Na - f0O_Na - fi5_Na)*C5_Na + C4_Na*f04_Na +
    I5_Na*bi5_Na + O_Na*b0O_Na)*allow_change;
  values[5] = ((-b0O_Na - fin_Na - fip_Na)*O_Na + B_Na*bip_Na + C5_Na*f0O_Na
    + I6_Na*bin_Na)*allow_change;
  values[6] = (O_Na*fip_Na - B_Na*bip_Na)*allow_change;
  values[7] = ((-bi1_Na - f11_Na)*I1_Na + C1_Na*fi1_Na +
    I2_Na*b11_Na)*allow_change;
  values[8] = ((-b11_Na - bi2_Na - f12_Na)*I2_Na + C2_Na*fi2_Na +
    I1_Na*f11_Na + I3_Na*b12_Na)*allow_change;
  values[9] = ((-b12_Na - bi3_Na - f13_Na)*I3_Na + C3_Na*fi3_Na +
    I2_Na*f12_Na + I4_Na*b13_Na)*allow_change;
  values[10] = ((-b13_Na - bi4_Na - f14_Na)*I4_Na + C4_Na*fi4_Na +
    I3_Na*f13_Na + I5_Na*b14_Na)*allow_change;
  values[11] = ((-b14_Na - bi5_Na - f1n_Na)*I5_Na + C5_Na*fi5_Na +
    I4_Na*f14_Na + I6_Na*b1n_Na)*allow_change;
  const double I_Na = g_Na*(-E_Na + v)*O_Na;

  // Expressions for the Kv1.1 component
  const double q10_Kv11 = 2.70000000000000;
  const double ca_Kv11 = 0.128890000000000;
  const double cva_Kv11 = 45.;
  const double cka_Kv11 = -33.9087700000000;
  const double cb_Kv11 = 0.128890000000000;
  const double cvb_Kv11 = 45.;
  const double ckb_Kv11 = 12.4210100000000;
  const double qt_Kv11 = std::pow(q10_Kv11, -2.2 + 0.1*Temp);
  const double alphan_Kv11 = ca_Kv11*std::exp((-cva_Kv11 - v)/cka_Kv11);
  const double betan_Kv11 = cb_Kv11*std::exp((-cvb_Kv11 - v)/ckb_Kv11);
  const double ninf_Kv11 = alphan_Kv11/(alphan_Kv11 + betan_Kv11);
  const double taun_Kv11 = 1./((alphan_Kv11 + betan_Kv11)*qt_Kv11);
  values[12] = (-n_Kv11 + ninf_Kv11)*allow_change/taun_Kv11;
  const double I_Kv11 = g_Kv11*std::pow(n_Kv11, 4.)*(-E_K + v);

  // Expressions for the Kv1.5 component
  const double q10_Kv15 = std::pow(2.2, -3.7 + 0.1*Temp);
  const double am_Kv15 =
    0.65*q10_Kv15/(1.66275285675459*std::exp(-0.0169491525423729*v) +
    0.308365167896581*std::exp(-0.117647058823529*v));
  const double bm_Kv15 = 0.65*q10_Kv15/(2.5 +
    124.403387639703*std::exp(0.0588235294117647*v));
  const double mtau_Kv15 = 1./(3.*(am_Kv15 + bm_Kv15));
  const double minf_Kv15 = 1.0/(1.0 +
    0.0425851362887876*std::exp(-0.104166666666667*v));
  const double an_Kv15 = 0.001*q10_Kv15/(2.4 +
    3.43809189356453*std::exp(-0.0128205128205128*v));
  const double bn_Kv15 = 2.75364493497472e-8*std::exp(0.0625*v)*q10_Kv15;
  const double ntau_Kv15 = 0.333333333333333/(an_Kv15 + bn_Kv15);
  const double ninf_Kv15 = 0.25 + 1.0/(1.35 +
    1.64872127070013*std::exp(0.0714285714285714*v));
  const double uinf_Kv15 = 0.1 + 1.0/(1.1 +
    1.64872127070013*std::exp(0.0714285714285714*v));
  const double utau_Kv15 = 6800.00000000000;
  values[13] = (-m_Kv15 + minf_Kv15)*allow_change/mtau_Kv15;
  values[14] = (-n_Kv15 + ninf_Kv15)*allow_change/ntau_Kv15;
  values[15] = (-u_Kv15 + uinf_Kv15)*allow_change/utau_Kv15;
  const double I_Kv15 = g_Kv15*(m_Kv15*m_Kv15*m_Kv15)*(0.1 + 1.0/(1.0 +
    3.17036319488807*std::exp(-0.0769230769230769*v)))*(-E_K +
    v)*n_Kv15*u_Kv15;

  // Expressions for the Kv3.3 component
  const double q10_Kv33 = 2.70000000000000;
  const double qt_Kv33 = std::pow(q10_Kv33, -2.2 + 0.1*Temp);
  const double ca_Kv33 = 0.220000000000000;
  const double cva_Kv33 = 16.;
  const double cka_Kv33 = -26.5000000000000;
  const double cb_Kv33 = 0.220000000000000;
  const double cvb_Kv33 = 16.;
  const double ckb_Kv33 = 26.5000000000000;
  const double alpha_Kv33 = ca_Kv33*std::exp((-cva_Kv33 - v)/cka_Kv33)*qt_Kv33;
  const double beta_Kv33 = cb_Kv33*std::exp((-cvb_Kv33 - v)/ckb_Kv33)*qt_Kv33;
  values[16] = ((1. - n_Kv33)*alpha_Kv33 - beta_Kv33*n_Kv33)*allow_change;
  const double I_Kv33 = g_Kv33*std::pow(n_Kv33, 4.)*(-E_K + v);

  // Expressions for the Kv3.4 component
  const double q10_Kv34 = 3.;
  const double qt_Kv34 = std::pow(q10_Kv34, -3.7 + 0.1*Temp);
  const double mivh_Kv34 = -24.;
  const double mik_Kv34 = 15.4000000000000;
  const double mty0_Kv34 = 0.000128510000000000;
  const double mtvh1_Kv34 = 100.700000000000;
  const double mtk1_Kv34 = 12.9000000000000;
  const double mtvh2_Kv34 = -56.0000000000000;
  const double mtk2_Kv34 = -23.1000000000000;
  const double hiy0_Kv34 = 0.310000000000000;
  const double hiA_Kv34 = 0.690000000000000;
  const double hivh_Kv34 = -5.80200000000000;
  const double hik_Kv34 = 11.2000000000000;
  const double minf_Kv34 = 1.0/(1.0 + std::exp((-11.0 + mivh_Kv34 -
    v)/mik_Kv34));
  const double mtau_Kv34 = 1000.*(11. + v < -35. ? 0.000102675 +
    0.0220402902836158*std::exp(0.0353481795687522*v) : mty0_Kv34 +
    1.0/(std::exp((11.0 + mtvh1_Kv34 + v)/mtk1_Kv34) + std::exp((11.0 +
    mtvh2_Kv34 + v)/mtk2_Kv34)))/qt_Kv34;
  const double hinf_Kv34 = hiy0_Kv34 + hiA_Kv34/(1. + std::exp((11.0 -
    hivh_Kv34 + v)/hik_Kv34));
  const double htau_Kv34 = 1000.*(11. + v > 0. ? 0.0012 +
    0.000487682413465536*std::exp(-0.141*v) : 1.2202e-5 +
    0.012*std::exp(-((1.35685483870968 +
    0.0201612903225806*v)*(1.35685483870968 +
    0.0201612903225806*v))))/qt_Kv34;
  values[17] = (-m_Kv34 + minf_Kv34)*allow_change/mtau_Kv34;
  values[18] = (-h_Kv34 + hinf_Kv34)*allow_change/htau_Kv34;
  const double I_Kv34 = g_Kv34*(m_Kv34*m_Kv34*m_Kv34)*(-E_K + v)*h_Kv34;

  // Expressions for the Kv4.3 component
  const double Q10_Kv43 = std::pow(3., -2.55 + 0.1*Temp);
  const double Aalpha_a_Kv43 = 0.814700000000000;
  const double Kalpha_a_Kv43 = -23.3270800000000;
  const double V0alpha_a_Kv43 = -9.17203000000000;
  const double Abeta_a_Kv43 = 0.165500000000000;
  const double Kbeta_a_Kv43 = 19.4717500000000;
  const double V0beta_a_Kv43 = -18.2791400000000;
  const double Aalpha_b_Kv43 = 0.0368000000000000;
  const double Kalpha_b_Kv43 = 12.8433000000000;
  const double V0alpha_b_Kv43 = -111.332090000000;
  const double Abeta_b_Kv43 = 0.0345000000000000;
  const double Kbeta_b_Kv43 = -8.90123000000000;
  const double V0beta_b_Kv43 = -49.9537000000000;
  const double alpha_a_Kv43 = Aalpha_a_Kv43*Q10_Kv43/(1. +
    std::exp((-V0alpha_a_Kv43 + v)/Kalpha_a_Kv43));
  const double beta_a_Kv43 = Abeta_a_Kv43*Q10_Kv43*std::exp(-(-V0beta_a_Kv43 +
    v)/Kbeta_a_Kv43);
  const double alpha_b_Kv43 = Aalpha_b_Kv43*Q10_Kv43/(1. +
    std::exp((-V0alpha_b_Kv43 + v)/Kalpha_b_Kv43));
  const double beta_b_Kv43 = Abeta_b_Kv43*Q10_Kv43/(1. +
    std::exp((-V0beta_b_Kv43 + v)/Kbeta_b_Kv43));
  const double a_inf_Kv43 = alpha_a_Kv43/(alpha_a_Kv43 + beta_a_Kv43);
  const double tau_a_Kv43 = 1.0/(alpha_a_Kv43 + beta_a_Kv43);
  const double b_inf_Kv43 = alpha_b_Kv43/(alpha_b_Kv43 + beta_b_Kv43);
  const double tau_b_Kv43 = 1.0/(alpha_b_Kv43 + beta_b_Kv43);
  values[19] = (-a_Kv43 + a_inf_Kv43)*allow_change/tau_a_Kv43;
  values[20] = (-b_Kv43 + b_inf_Kv43)*allow_change/tau_b_Kv43;
  const double I_Kv43 = g_Kv43*(a_Kv43*a_Kv43*a_Kv43)*(-E_K + v)*b_Kv43;

  // Expressions for the Kir2.x component
  const double Q10_Kir2x = std::pow(3., -2.0 + 0.1*Temp);
  const double Aalpha_d_Kir2x = 0.132890000000000;
  const double Kalpha_d_Kir2x = -24.3902000000000;
  const double V0alpha_d_Kir2x = -83.9400000000000;
  const double Abeta_d_Kir2x = 0.169940000000000;
  const double Kbeta_d_Kir2x = 35.7140000000000;
  const double V0beta_d_Kir2x = -83.9400000000000;
  const double a_d_Kir2x =
    Aalpha_d_Kir2x*Q10_Kir2x*std::exp((-V0alpha_d_Kir2x + v)/Kalpha_d_Kir2x);
  const double b_d_Kir2x = Abeta_d_Kir2x*Q10_Kir2x*std::exp((-V0beta_d_Kir2x
    + v)/Kbeta_d_Kir2x);
  const double tau_d_Kir2x = 1.0/(a_d_Kir2x + b_d_Kir2x);
  const double d_inf_Kir2x = a_d_Kir2x/(a_d_Kir2x + b_d_Kir2x);
  values[21] = (-d_Kir2x + d_inf_Kir2x)*allow_change/tau_d_Kir2x;
  const double I_Kir2x = g_Kir2x*(-E_K + v)*d_Kir2x;

  // Expressions for the Kca1.1 component
  const double q10_Kca11 = 3.;
  const double qt_Kca11 = std::pow(q10_Kca11, -2.3 + 0.1*Temp);
  const double Qo_Kca11 = 0.730000000000000;
  const double Qc_Kca11 = -0.670000000000000;
  const double k1_Kca11 = 1000.00000000000;
  const double onoffrate_Kca11 = 1.;
  const double Kc_Kca11 = 0.0110000000000000;
  const double Ko_Kca11 = 0.00110000000000000;
  const double pf0_Kca11 = 0.00239000000000000;
  const double pf1_Kca11 = 0.00700000000000000;
  const double pf2_Kca11 = 0.0400000000000000;
  const double pf3_Kca11 = 0.295000000000000;
  const double pf4_Kca11 = 0.557000000000000;
  const double pb0_Kca11 = 3.93600000000000;
  const double pb1_Kca11 = 1.15200000000000;
  const double pb2_Kca11 = 0.659000000000000;
  const double pb3_Kca11 = 0.486000000000000;
  const double pb4_Kca11 = 0.0920000000000000;
  const double c01_Kca11 = 4.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c12_Kca11 = 3.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c23_Kca11 = 2.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c34_Kca11 = cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o01_Kca11 = 4.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o12_Kca11 = 3.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o23_Kca11 = 2.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o34_Kca11 = cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c10_Kca11 = Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c21_Kca11 = 2.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c32_Kca11 = 3.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c43_Kca11 = 4.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o10_Kca11 = Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o21_Kca11 = 2.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o32_Kca11 = 3.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o43_Kca11 = 4.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double alpha_Kca11 = std::exp(FARADAY*Qo_Kca11*v/(R*(273.15 + Temp)));
  const double beta_Kca11 = std::exp(FARADAY*Qc_Kca11*v/(R*(273.15 + Temp)));
  const double f0_Kca11 = pf0_Kca11*alpha_Kca11*qt_Kca11;
  const double f1_Kca11 = pf1_Kca11*alpha_Kca11*qt_Kca11;
  const double f2_Kca11 = pf2_Kca11*alpha_Kca11*qt_Kca11;
  const double f3_Kca11 = pf3_Kca11*alpha_Kca11*qt_Kca11;
  const double f4_Kca11 = pf4_Kca11*alpha_Kca11*qt_Kca11;
  const double b0_Kca11 = pb0_Kca11*beta_Kca11*qt_Kca11;
  const double b1_Kca11 = pb1_Kca11*beta_Kca11*qt_Kca11;
  const double b2_Kca11 = pb2_Kca11*beta_Kca11*qt_Kca11;
  const double b3_Kca11 = pb3_Kca11*beta_Kca11*qt_Kca11;
  const double b4_Kca11 = pb4_Kca11*beta_Kca11*qt_Kca11;
  const double C4_Kca11 = 1. - C0_Kca11 - C1_Kca11 - C2_Kca11 - C3_Kca11 -
    O0_Kca11 - O1_Kca11 - O2_Kca11 - O3_Kca11 - O4_Kca11;
  values[22] = ((-c01_Kca11 - f0_Kca11)*C0_Kca11 + C1_Kca11*c10_Kca11 +
    O0_Kca11*b0_Kca11)*allow_change;
  values[23] = ((-c10_Kca11 - c12_Kca11 - f1_Kca11)*C1_Kca11 +
    C0_Kca11*c01_Kca11 + C2_Kca11*c21_Kca11 + O1_Kca11*b1_Kca11)*allow_change;
  values[24] = ((-c21_Kca11 - c23_Kca11 - f2_Kca11)*C2_Kca11 +
    C1_Kca11*c12_Kca11 + C3_Kca11*c32_Kca11 + O2_Kca11*b2_Kca11)*allow_change;
  values[25] = ((-c32_Kca11 - c34_Kca11 - f3_Kca11)*C3_Kca11 +
    C2_Kca11*c23_Kca11 + C4_Kca11*c43_Kca11 + O3_Kca11*b3_Kca11)*allow_change;
  values[26] = ((-b0_Kca11 - o01_Kca11)*O0_Kca11 + C0_Kca11*f0_Kca11 +
    O1_Kca11*o10_Kca11)*allow_change;
  values[27] = ((-b1_Kca11 - o10_Kca11 - o12_Kca11)*O1_Kca11 +
    C1_Kca11*f1_Kca11 + O0_Kca11*o01_Kca11 + O2_Kca11*o21_Kca11)*allow_change;
  values[28] = ((-b2_Kca11 - o21_Kca11 - o23_Kca11)*O2_Kca11 +
    C2_Kca11*f2_Kca11 + O1_Kca11*o12_Kca11 + O3_Kca11*o32_Kca11)*allow_change;
  values[29] = ((-b3_Kca11 - o32_Kca11 - o34_Kca11)*O3_Kca11 +
    C3_Kca11*f3_Kca11 + O2_Kca11*o23_Kca11 + O4_Kca11*o43_Kca11)*allow_change;
  values[30] = ((-b4_Kca11 - o43_Kca11)*O4_Kca11 + C4_Kca11*f4_Kca11 +
    O3_Kca11*o34_Kca11)*allow_change;
  const double I_Kca11 = g_Kca11*(-E_K + v)*(O0_Kca11 + O1_Kca11 + O2_Kca11 +
    O3_Kca11 + O4_Kca11);

  // Expressions for the Kca2.2 component
  const double Q10_Kca22 = 3.00000000000000;
  const double tcorr_Kca22 = std::pow(Q10_Kca22, -2.3 + 0.1*Temp);
  const double invc1_Kca22 = 0.0800000000000000;
  const double invc2_Kca22 = 0.0800000000000000;
  const double invc3_Kca22 = 0.200000000000000;
  const double invo1_Kca22 = 1.;
  const double invo2_Kca22 = 0.100000000000000;
  const double diro1_Kca22 = 0.160000000000000;
  const double diro2_Kca22 = 1.20000000000000;
  const double dirc2_Kca22 = 200.000000000000;
  const double dirc3_Kca22 = 160.000000000000;
  const double dirc4_Kca22 = 80.0000000000000;
  const double invc1_t_Kca22 = invc1_Kca22*tcorr_Kca22;
  const double invc2_t_Kca22 = invc2_Kca22*tcorr_Kca22;
  const double invc3_t_Kca22 = invc3_Kca22*tcorr_Kca22;
  const double invo1_t_Kca22 = invo1_Kca22*tcorr_Kca22;
  const double invo2_t_Kca22 = invo2_Kca22*tcorr_Kca22;
  const double diro1_t_Kca22 = diro1_Kca22*tcorr_Kca22;
  const double diro2_t_Kca22 = diro2_Kca22*tcorr_Kca22;
  const double dirc2_t_Kca22 = dirc2_Kca22*tcorr_Kca22;
  const double dirc3_t_Kca22 = dirc3_Kca22*tcorr_Kca22;
  const double dirc4_t_Kca22 = dirc4_Kca22*tcorr_Kca22;
  const double dirc2_t_ca_Kca22 = cai*dirc2_t_Kca22;
  const double dirc3_t_ca_Kca22 = cai*dirc3_t_Kca22;
  const double dirc4_t_ca_Kca22 = cai*dirc4_t_Kca22;
  const double c4_Kca22 = 1. - c1_Kca22 - c2_Kca22 - c3_Kca22 - o1_Kca22 -
    o2_Kca22;
  values[31] = (c2_Kca22*invc1_t_Kca22 -
    c1_Kca22*dirc2_t_ca_Kca22)*allow_change;
  values[32] = ((-dirc3_t_ca_Kca22 - invc1_t_Kca22)*c2_Kca22 +
    c1_Kca22*dirc2_t_ca_Kca22 + c3_Kca22*invc2_t_Kca22)*allow_change;
  values[33] = ((-dirc4_t_ca_Kca22 - diro1_t_Kca22 - invc2_t_Kca22)*c3_Kca22 +
    c2_Kca22*dirc3_t_ca_Kca22 + c4_Kca22*invc3_t_Kca22 +
    invo1_t_Kca22*o1_Kca22)*allow_change;
  values[34] = (c3_Kca22*diro1_t_Kca22 - invo1_t_Kca22*o1_Kca22)*allow_change;
  values[35] = (c4_Kca22*diro2_t_Kca22 - invo2_t_Kca22*o2_Kca22)*allow_change;
  const double I_Kca22 = g_Kca22*(-E_K + v)*(o1_Kca22 + o2_Kca22);

  // Expressions for the Kca3.1 component
  const double Yvdep_Kca31 = 13.364375107327*std::exp(0.037037037037037*v);
  const double Yconcdep_Kca31 = (cai < 0.01 ? (7.5 - 500.*cai)/(-1. +
    102586.491213197*std::exp(-769.230769230769*cai)) : 0.0545700593112901);
  const double Ybeta_Kca31 = 0.0500000000000000;
  const double Yalpha_Kca31 = Yconcdep_Kca31*Yvdep_Kca31;
  values[36] = ((1. - Y_Kca31)*Yalpha_Kca31 -
    Ybeta_Kca31*Y_Kca31)*allow_change;
  const double I_Kca31 = g_Kca31*(-E_K + v)*Y_Kca31;

  // Expressions for the Cav2.1 component
  const double q10_Cav21 = 3.;
  const double qt_Cav21 = std::pow(q10_Cav21, -2.3 + 0.1*Temp);
  const double vhalfm_Cav21 = -29.4580000000000;
  const double cvm_Cav21 = 8.42900000000000;
  const double vshift_Cav21 = 0.;
  const double minf_Cav21 = 1.0/(1.0 + std::exp((vhalfm_Cav21 + vshift_Cav21 -
    v)/cvm_Cav21));
  const double taum_Cav21 = 1.0*(v >= -40. ? 0.2702 +
    1.1622*std::exp(0.00609050490285645*(26.798 + v)*(-26.798 - v)) :
    0.6923*std::exp(0.00091796007240869*v))/qt_Cav21;
  values[37] = (-m_Cav21 + minf_Cav21)*allow_change/taum_Cav21;
  const double z_Ca = 2.;
  const double zeta_Ca = FARADAY*z_Ca*v/(R*(273.19 + Temp));
  const double ghk = (std::fabs(1. - std::exp(-zeta_Ca)) < 1.0e-6 ?
    1.0e-6*F*z_Ca*(1. + zeta_Ca/2.)*(cai - cao*std::exp(-zeta_Ca)) :
    1.0e-6*F*z_Ca*(cai - cao*std::exp(-zeta_Ca))*zeta_Ca/(1. -
    std::exp(-zeta_Ca)));
  const double I_Cav21 = 1000.0*g_Cav21*(m_Cav21*m_Cav21*m_Cav21)*ghk;

  // Expressions for the Cav3.1 component
  const double q10_Cav31 = 3.;
  const double qt_Cav31 = std::pow(q10_Cav31, -3.7 + 0.1*Temp);
  const double v0_m_inf_Cav31 = -52.0000000000000;
  const double v0_h_inf_Cav31 = -72.0000000000000;
  const double k_m_inf_Cav31 = -5.00000000000000;
  const double k_h_inf_Cav31 = 7.00000000000000;
  const double C_tau_m_Cav31 = 1.00000000000000;
  const double A_tau_m_Cav31 = 1.00000000000000;
  const double v0_tau_m1_Cav31 = -40.0000000000000;
  const double v0_tau_m2_Cav31 = -102.000000000000;
  const double k_tau_m1_Cav31 = 9.00000000000000;
  const double k_tau_m2_Cav31 = -18.0000000000000;
  const double C_tau_h_Cav31 = 15.0000000000000;
  const double A_tau_h_Cav31 = 1.00000000000000;
  const double v0_tau_h1_Cav31 = -32.0000000000000;
  const double k_tau_h1_Cav31 = 7.00000000000000;
  const double minf_Cav31 = 1.0/(1. + std::exp((-v0_m_inf_Cav31 +
    v)/k_m_inf_Cav31));
  const double hinf_Cav31 = 1.0/(1. + std::exp((-v0_h_inf_Cav31 +
    v)/k_h_inf_Cav31));
  const double taum_Cav31 = (v <= -90. ? 1. : (C_tau_m_Cav31 +
    A_tau_m_Cav31/(std::exp((-v0_tau_m1_Cav31 + v)/k_tau_m1_Cav31) +
    std::exp((-v0_tau_m2_Cav31 + v)/k_tau_m2_Cav31)))/qt_Cav31);
  const double tauh_Cav31 = (C_tau_h_Cav31 +
    A_tau_h_Cav31*std::exp(-(-v0_tau_h1_Cav31 + v)/k_tau_h1_Cav31))/qt_Cav31;
  values[38] = (-m_Cav31 + minf_Cav31)*allow_change/taum_Cav31;
  values[39] = (-h_Cav31 + hinf_Cav31)*allow_change/tauh_Cav31;
  const double I_Cav31 = 1000.0*g_Cav31*(m_Cav31*m_Cav31)*ghk*h_Cav31;

  // Expressions for the Cav3.2 component
  const double phi_m_Cav32 = 6.89864830730607;
  const double phi_h_Cav32 = 3.73719281884655;
  const double m_inf_Cav32 = 1.0/(1. +
    0.000607957605998435*std::exp(-0.135135135135135*v));
  const double h_inf_Cav32 = 1.0/(1. +
    148461.062050627*std::exp(0.139275766016713*v));
  const double tau_m_Cav32 = (1.9 +
    1.0/(22.404093719072*std::exp(0.0840336134453781*v) +
    0.00189854653589182*std::exp(-v/21.)))/phi_m_Cav32;
  const double tau_h_Cav32 = 13.7 + (1942.0 +
    55178666.3027343*std::exp(0.108695652173913*v))/(phi_h_Cav32*(1. +
    30321871936.1927*std::exp(0.27027027027027*v)));
  values[40] = (-m_Cav32 + m_inf_Cav32)*allow_change/tau_m_Cav32;
  values[41] = (-h_Cav32 + h_inf_Cav32)*allow_change/tau_h_Cav32;
  const double E_Ca = R*(273.15 + Temp)*std::log(cao/cai)/(2.*FARADAY);
  const double I_Cav32 = g_Cav32*(m_Cav32*m_Cav32)*(-E_Ca + v)*h_Cav32;

  // Expressions for the Cav3.3 component
  const double q10_Cav33 = 2.30000000000000;
  const double qt_Cav33 = std::pow(q10_Cav33, -2.8 + 0.1*Temp);
  const double gCav3_3bar = 1.00000000000000e-5;
  const double vhalfn_Cav33 = -41.5000000000000;
  const double vhalfl_Cav33 = -69.8000000000000;
  const double kn_Cav33 = 6.20000000000000;
  const double kl_Cav33 = -6.10000000000000;
  const double n_inf_Cav33 = 1.0/(1.0 + std::exp((vhalfn_Cav33 - v)/kn_Cav33));
  const double l_inf_Cav33 = 1.0/(1.0 + std::exp((vhalfl_Cav33 - v)/kl_Cav33));
  const double tau_n_Cav33 = (v > -60. ? 7.2 +
    0.02*std::exp(-0.0680272108843537*v) : 79.5 +
    2.0*std::exp(-0.10752688172043*v))/qt_Cav33;
  const double tau_l_Cav33 = (v > -60. ? 0.875*std::exp(120./41. + v/41.) :
    260.0)/qt_Cav33;
  const double w_Cav33 = FARADAY*z_Ca*v/(R*(273.14 + Temp));
  const double ghk_Cav33 = -FARADAY*z_Ca*(cao -
    cai*std::exp(w_Cav33))*w_Cav33/(-1.0 + std::exp(w_Cav33));
  values[42] = (-n_Cav33 + n_inf_Cav33)*allow_change/tau_n_Cav33;
  values[43] = (-l_Cav33 + l_inf_Cav33)*allow_change/tau_l_Cav33;
  const double I_Cav33 =
    gCav3_3bar*g_Cav33*(n_Cav33*n_Cav33*n_Cav33)*ghk_Cav33*l_Cav33;

  // Expressions for the HCN1 component
  const double q10_HCN1 = 3.00000000000000;
  const double qt_HCN1 = std::pow(q10_HCN1, -3.7 + 0.1*Temp);
  const double ratetau_HCN1 = 1.00000000000000;
  const double ljp_HCN1 = 9.30000000000000;
  const double v_inf_half_noljp_HCN1 = -90.3000000000000;
  const double v_inf_k_HCN1 = 9.67000000000000;
  const double v_tau_const_HCN1 = 0.00180000000000000;
  const double v_tau_half1_noljp_HCN1 = -68.0000000000000;
  const double v_tau_half2_noljp_HCN1 = -68.0000000000000;
  const double v_tau_k1_HCN1 = -22.0000000000000;
  const double v_tau_k2_HCN1 = 7.14000000000000;
  const double v_inf_half_HCN1 = v_inf_half_noljp_HCN1 - ljp_HCN1;
  const double v_tau_half1_HCN1 = v_tau_half1_noljp_HCN1 - ljp_HCN1;
  const double v_tau_half2_HCN1 = v_tau_half2_noljp_HCN1 - ljp_HCN1;
  const double hinf_HCN1 = 1.0/(1. + std::exp((-v_inf_half_HCN1 +
    v)/v_inf_k_HCN1));
  const double tauh_HCN1 =
    ratetau_HCN1/(v_tau_const_HCN1*(std::exp((-v_tau_half1_HCN1 +
    v)/v_tau_k1_HCN1) + std::exp((-v_tau_half2_HCN1 +
    v)/v_tau_k2_HCN1))*qt_HCN1);
  values[44] = (-h_HCN1 + hinf_HCN1)*allow_change/tauh_HCN1;
  const double I_HCN1 = g_HCN1*(-E_h + v)*h_HCN1;

  // Expressions for the leak component
  const double I_leak = g_leak*(-E_leak + v);

  // Expressions for the Membrane component
  const double I_tot = I_Cav21 + I_Cav31 + I_Cav32 + I_Cav33 + I_HCN1 +
    I_Kca11 + I_Kca22 + I_Kca31 + I_Kir2x + I_Kv11 + I_Kv15 + I_Kv33 + I_Kv34 +
    I_Kv43 + I_Na + I_leak;
  values[45] = -1.0*I_tot*allow_change/Cm;
}

// Compute a forward step using the rush larsen algorithm to the Masoli2015 ODE
void forward_rush_larsen(double* states, const double t, const double dt,
  const double* parameters)
{

  // Assign states
  const double C1_Na = states[0];
  const double C2_Na = states[1];
  const double C3_Na = states[2];
  const double C4_Na = states[3];
  const double C5_Na = states[4];
  const double O_Na = states[5];
  const double B_Na = states[6];
  const double I1_Na = states[7];
  const double I2_Na = states[8];
  const double I3_Na = states[9];
  const double I4_Na = states[10];
  const double I5_Na = states[11];
  const double n_Kv11 = states[12];
  const double m_Kv15 = states[13];
  const double n_Kv15 = states[14];
  const double u_Kv15 = states[15];
  const double n_Kv33 = states[16];
  const double m_Kv34 = states[17];
  const double h_Kv34 = states[18];
  const double a_Kv43 = states[19];
  const double b_Kv43 = states[20];
  const double d_Kir2x = states[21];
  const double C0_Kca11 = states[22];
  const double C1_Kca11 = states[23];
  const double C2_Kca11 = states[24];
  const double C3_Kca11 = states[25];
  const double O0_Kca11 = states[26];
  const double O1_Kca11 = states[27];
  const double O2_Kca11 = states[28];
  const double O3_Kca11 = states[29];
  const double O4_Kca11 = states[30];
  const double c1_Kca22 = states[31];
  const double c2_Kca22 = states[32];
  const double c3_Kca22 = states[33];
  const double o1_Kca22 = states[34];
  const double o2_Kca22 = states[35];
  const double Y_Kca31 = states[36];
  const double m_Cav21 = states[37];
  const double m_Cav31 = states[38];
  const double h_Cav31 = states[39];
  const double m_Cav32 = states[40];
  const double h_Cav32 = states[41];
  const double n_Cav33 = states[42];
  const double l_Cav33 = states[43];
  const double h_HCN1 = states[44];
  const double v = states[45];

  // Assign parameters
  const double Cm = parameters[0];
  const double F = parameters[1];
  const double FARADAY = parameters[2];
  const double R = parameters[3];
  const double Tdiff = parameters[4];
  const double Temp = parameters[5];
  const double cai = parameters[6];
  const double cao = parameters[7];
  const double E_Na = parameters[8];
  const double g_Na = parameters[9];
  const double E_K = parameters[10];
  const double g_Kv11 = parameters[11];
  const double g_Kv15 = parameters[12];
  const double g_Kv33 = parameters[13];
  const double g_Kv34 = parameters[14];
  const double g_Kv43 = parameters[15];
  const double g_Kir2x = parameters[16];
  const double g_Kca11 = parameters[17];
  const double g_Kca22 = parameters[18];
  const double g_Kca31 = parameters[19];
  const double g_Cav21 = parameters[20];
  const double g_Cav31 = parameters[21];
  const double g_Cav32 = parameters[22];
  const double g_Cav33 = parameters[23];
  const double E_h = parameters[24];
  const double g_HCN1 = parameters[25];
  const double E_leak = parameters[26];
  const double g_leak = parameters[27];

  // Expressions for the Nav1.6 component
  const double allow_change = (t >= Tdiff ? 1. : 0.);
  const double q10_Na = 3.;
  const double qt_Na = std::pow(q10_Na, -2.2 + 0.1*Temp);
  const double Con_Na = 0.00500000000000000;
  const double Coff_Na = 0.500000000000000;
  const double Oon_Na = 0.750000000000000;
  const double Ooff_Na = 0.00500000000000000;
  const double alpha_Na = 150.;
  const double beta_Na = 3.;
  const double gamma_Na = 150.;
  const double delta_Na = 40.;
  const double epsilon_Na = 1.75000000000000;
  const double zeta_Na = 0.0300000000000000;
  const double x1_Na = 20.;
  const double x2_Na = -20.;
  const double x3_Na = 1000000000000.00;
  const double x4_Na = -1000000000000.00;
  const double x5_Na = 1000000000000.00;
  const double x6_Na = -25.;
  const double alfac_Na = std::pow(Oon_Na/Con_Na, 0.25);
  const double btfac_Na = std::pow(Ooff_Na/Coff_Na, 0.25);
  const double f01_Na = 4.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f02_Na = 3.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f03_Na = 2.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f04_Na = alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f0O_Na = gamma_Na*std::exp(v/x3_Na)*qt_Na;
  const double fip_Na = epsilon_Na*std::exp(v/x5_Na)*qt_Na;
  const double f11_Na = 4.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f12_Na = 3.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f13_Na = 2.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f14_Na = alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f1n_Na = gamma_Na*std::exp(v/x3_Na)*qt_Na;
  const double fi1_Na = Con_Na*qt_Na;
  const double fi2_Na = Con_Na*alfac_Na*qt_Na;
  const double fi3_Na = Con_Na*(alfac_Na*alfac_Na)*qt_Na;
  const double fi4_Na = Con_Na*(alfac_Na*alfac_Na*alfac_Na)*qt_Na;
  const double fi5_Na = Con_Na*std::pow(alfac_Na, 4.)*qt_Na;
  const double fin_Na = Oon_Na*qt_Na;
  const double b01_Na = beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b02_Na = 2.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b03_Na = 3.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b04_Na = 4.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b0O_Na = delta_Na*std::exp(v/x4_Na)*qt_Na;
  const double bip_Na = zeta_Na*std::exp(v/x6_Na)*qt_Na;
  const double b11_Na = beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b12_Na = 2.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b13_Na = 3.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b14_Na = 4.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b1n_Na = delta_Na*std::exp(v/x4_Na)*qt_Na;
  const double bi1_Na = Coff_Na*qt_Na;
  const double bi2_Na = Coff_Na*btfac_Na*qt_Na;
  const double bi3_Na = Coff_Na*(btfac_Na*btfac_Na)*qt_Na;
  const double bi4_Na = Coff_Na*(btfac_Na*btfac_Na*btfac_Na)*qt_Na;
  const double bi5_Na = Coff_Na*std::pow(btfac_Na, 4.)*qt_Na;
  const double bin_Na = Ooff_Na*qt_Na;
  const double I6_Na = 1. - B_Na - C1_Na - C2_Na - C3_Na - C4_Na - C5_Na -
    I1_Na - I2_Na - I3_Na - I4_Na - I5_Na - O_Na;
  const double dC1_Na_dt = ((-f01_Na - fi1_Na)*C1_Na + C2_Na*b01_Na +
    I1_Na*bi1_Na)*allow_change;
  const double dC1_Na_dt_linearized = (-f01_Na - fi1_Na)*allow_change;
  states[0] = C1_Na + (std::fabs(dC1_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dC1_Na_dt_linearized))*dC1_Na_dt/dC1_Na_dt_linearized :
    dt*dC1_Na_dt);
  const double dC2_Na_dt = ((-b01_Na - f02_Na - fi2_Na)*C2_Na + C1_Na*f01_Na
    + C3_Na*b02_Na + I2_Na*bi2_Na)*allow_change;
  const double dC2_Na_dt_linearized = (-b01_Na - f02_Na - fi2_Na)*allow_change;
  states[1] = C2_Na + (std::fabs(dC2_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dC2_Na_dt_linearized))*dC2_Na_dt/dC2_Na_dt_linearized :
    dt*dC2_Na_dt);
  const double dC3_Na_dt = ((-b02_Na - f03_Na - fi3_Na)*C3_Na + C2_Na*f02_Na
    + C4_Na*b03_Na + I3_Na*bi3_Na)*allow_change;
  const double dC3_Na_dt_linearized = (-b02_Na - f03_Na - fi3_Na)*allow_change;
  states[2] = C3_Na + (std::fabs(dC3_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dC3_Na_dt_linearized))*dC3_Na_dt/dC3_Na_dt_linearized :
    dt*dC3_Na_dt);
  const double dC4_Na_dt = ((-b03_Na - f04_Na - fi4_Na)*C4_Na + C3_Na*f03_Na
    + C5_Na*b04_Na + I4_Na*bi4_Na)*allow_change;
  const double dC4_Na_dt_linearized = (-b03_Na - f04_Na - fi4_Na)*allow_change;
  states[3] = C4_Na + (std::fabs(dC4_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dC4_Na_dt_linearized))*dC4_Na_dt/dC4_Na_dt_linearized :
    dt*dC4_Na_dt);
  const double dC5_Na_dt = ((-b04_Na - f0O_Na - fi5_Na)*C5_Na + C4_Na*f04_Na
    + I5_Na*bi5_Na + O_Na*b0O_Na)*allow_change;
  const double dC5_Na_dt_linearized = (-b04_Na - f0O_Na - fi5_Na)*allow_change;
  states[4] = C5_Na + (std::fabs(dC5_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dC5_Na_dt_linearized))*dC5_Na_dt/dC5_Na_dt_linearized :
    dt*dC5_Na_dt);
  const double dO_Na_dt = ((-b0O_Na - fin_Na - fip_Na)*O_Na + B_Na*bip_Na +
    C5_Na*f0O_Na + I6_Na*bin_Na)*allow_change;
  const double dO_Na_dt_linearized = (-b0O_Na - bin_Na - fin_Na -
    fip_Na)*allow_change;
  states[5] = O_Na + (std::fabs(dO_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dO_Na_dt_linearized))*dO_Na_dt/dO_Na_dt_linearized :
    dt*dO_Na_dt);
  const double dB_Na_dt = (O_Na*fip_Na - B_Na*bip_Na)*allow_change;
  const double dB_Na_dt_linearized = -allow_change*bip_Na;
  states[6] = B_Na + (std::fabs(dB_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dB_Na_dt_linearized))*dB_Na_dt/dB_Na_dt_linearized :
    dt*dB_Na_dt);
  const double dI1_Na_dt = ((-bi1_Na - f11_Na)*I1_Na + C1_Na*fi1_Na +
    I2_Na*b11_Na)*allow_change;
  const double dI1_Na_dt_linearized = (-bi1_Na - f11_Na)*allow_change;
  states[7] = I1_Na + (std::fabs(dI1_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dI1_Na_dt_linearized))*dI1_Na_dt/dI1_Na_dt_linearized :
    dt*dI1_Na_dt);
  const double dI2_Na_dt = ((-b11_Na - bi2_Na - f12_Na)*I2_Na + C2_Na*fi2_Na
    + I1_Na*f11_Na + I3_Na*b12_Na)*allow_change;
  const double dI2_Na_dt_linearized = (-b11_Na - bi2_Na - f12_Na)*allow_change;
  states[8] = I2_Na + (std::fabs(dI2_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dI2_Na_dt_linearized))*dI2_Na_dt/dI2_Na_dt_linearized :
    dt*dI2_Na_dt);
  const double dI3_Na_dt = ((-b12_Na - bi3_Na - f13_Na)*I3_Na + C3_Na*fi3_Na
    + I2_Na*f12_Na + I4_Na*b13_Na)*allow_change;
  const double dI3_Na_dt_linearized = (-b12_Na - bi3_Na - f13_Na)*allow_change;
  states[9] = I3_Na + (std::fabs(dI3_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dI3_Na_dt_linearized))*dI3_Na_dt/dI3_Na_dt_linearized :
    dt*dI3_Na_dt);
  const double dI4_Na_dt = ((-b13_Na - bi4_Na - f14_Na)*I4_Na + C4_Na*fi4_Na
    + I3_Na*f13_Na + I5_Na*b14_Na)*allow_change;
  const double dI4_Na_dt_linearized = (-b13_Na - bi4_Na - f14_Na)*allow_change;
  states[10] = I4_Na + (std::fabs(dI4_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dI4_Na_dt_linearized))*dI4_Na_dt/dI4_Na_dt_linearized :
    dt*dI4_Na_dt);
  const double dI5_Na_dt = ((-b14_Na - bi5_Na - f1n_Na)*I5_Na + C5_Na*fi5_Na
    + I4_Na*f14_Na + I6_Na*b1n_Na)*allow_change;
  const double dI5_Na_dt_linearized = (-b14_Na - b1n_Na - bi5_Na -
    f1n_Na)*allow_change;
  states[11] = I5_Na + (std::fabs(dI5_Na_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dI5_Na_dt_linearized))*dI5_Na_dt/dI5_Na_dt_linearized :
    dt*dI5_Na_dt);
  const double I_Na = g_Na*(-E_Na + v)*O_Na;

  // Expressions for the Kv1.1 component
  const double q10_Kv11 = 2.70000000000000;
  const double ca_Kv11 = 0.128890000000000;
  const double cva_Kv11 = 45.;
  const double cka_Kv11 = -33.9087700000000;
  const double cb_Kv11 = 0.128890000000000;
  const double cvb_Kv11 = 45.;
  const double ckb_Kv11 = 12.4210100000000;
  const double qt_Kv11 = std::pow(q10_Kv11, -2.2 + 0.1*Temp);
  const double alphan_Kv11 = ca_Kv11*std::exp((-cva_Kv11 - v)/cka_Kv11);
  const double betan_Kv11 = cb_Kv11*std::exp((-cvb_Kv11 - v)/ckb_Kv11);
  const double ninf_Kv11 = alphan_Kv11/(alphan_Kv11 + betan_Kv11);
  const double taun_Kv11 = 1./((alphan_Kv11 + betan_Kv11)*qt_Kv11);
  const double dn_Kv11_dt = (-n_Kv11 + ninf_Kv11)*allow_change/taun_Kv11;
  const double dn_Kv11_dt_linearized = -allow_change/taun_Kv11;
  states[12] = (std::fabs(dn_Kv11_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dn_Kv11_dt_linearized))*dn_Kv11_dt/dn_Kv11_dt_linearized :
    dt*dn_Kv11_dt) + n_Kv11;
  const double I_Kv11 = g_Kv11*std::pow(n_Kv11, 4.)*(-E_K + v);

  // Expressions for the Kv1.5 component
  const double q10_Kv15 = std::pow(2.2, -3.7 + 0.1*Temp);
  const double am_Kv15 =
    0.65*q10_Kv15/(1.66275285675459*std::exp(-0.0169491525423729*v) +
    0.308365167896581*std::exp(-0.117647058823529*v));
  const double bm_Kv15 = 0.65*q10_Kv15/(2.5 +
    124.403387639703*std::exp(0.0588235294117647*v));
  const double mtau_Kv15 = 1./(3.*(am_Kv15 + bm_Kv15));
  const double minf_Kv15 = 1.0/(1.0 +
    0.0425851362887876*std::exp(-0.104166666666667*v));
  const double an_Kv15 = 0.001*q10_Kv15/(2.4 +
    3.43809189356453*std::exp(-0.0128205128205128*v));
  const double bn_Kv15 = 2.75364493497472e-8*std::exp(0.0625*v)*q10_Kv15;
  const double ntau_Kv15 = 0.333333333333333/(an_Kv15 + bn_Kv15);
  const double ninf_Kv15 = 0.25 + 1.0/(1.35 +
    1.64872127070013*std::exp(0.0714285714285714*v));
  const double uinf_Kv15 = 0.1 + 1.0/(1.1 +
    1.64872127070013*std::exp(0.0714285714285714*v));
  const double utau_Kv15 = 6800.00000000000;
  const double dm_Kv15_dt = (-m_Kv15 + minf_Kv15)*allow_change/mtau_Kv15;
  const double dm_Kv15_dt_linearized = -allow_change/mtau_Kv15;
  states[13] = (std::fabs(dm_Kv15_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dm_Kv15_dt_linearized))*dm_Kv15_dt/dm_Kv15_dt_linearized :
    dt*dm_Kv15_dt) + m_Kv15;
  const double dn_Kv15_dt = (-n_Kv15 + ninf_Kv15)*allow_change/ntau_Kv15;
  const double dn_Kv15_dt_linearized = -allow_change/ntau_Kv15;
  states[14] = (std::fabs(dn_Kv15_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dn_Kv15_dt_linearized))*dn_Kv15_dt/dn_Kv15_dt_linearized :
    dt*dn_Kv15_dt) + n_Kv15;
  const double du_Kv15_dt = (-u_Kv15 + uinf_Kv15)*allow_change/utau_Kv15;
  const double du_Kv15_dt_linearized = -allow_change/utau_Kv15;
  states[15] = (std::fabs(du_Kv15_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*du_Kv15_dt_linearized))*du_Kv15_dt/du_Kv15_dt_linearized :
    dt*du_Kv15_dt) + u_Kv15;
  const double I_Kv15 = g_Kv15*(m_Kv15*m_Kv15*m_Kv15)*(0.1 + 1.0/(1.0 +
    3.17036319488807*std::exp(-0.0769230769230769*v)))*(-E_K +
    v)*n_Kv15*u_Kv15;

  // Expressions for the Kv3.3 component
  const double q10_Kv33 = 2.70000000000000;
  const double qt_Kv33 = std::pow(q10_Kv33, -2.2 + 0.1*Temp);
  const double ca_Kv33 = 0.220000000000000;
  const double cva_Kv33 = 16.;
  const double cka_Kv33 = -26.5000000000000;
  const double cb_Kv33 = 0.220000000000000;
  const double cvb_Kv33 = 16.;
  const double ckb_Kv33 = 26.5000000000000;
  const double alpha_Kv33 = ca_Kv33*std::exp((-cva_Kv33 - v)/cka_Kv33)*qt_Kv33;
  const double beta_Kv33 = cb_Kv33*std::exp((-cvb_Kv33 - v)/ckb_Kv33)*qt_Kv33;
  const double dn_Kv33_dt = ((1. - n_Kv33)*alpha_Kv33 -
    beta_Kv33*n_Kv33)*allow_change;
  const double dn_Kv33_dt_linearized = (-alpha_Kv33 - beta_Kv33)*allow_change;
  states[16] = (std::fabs(dn_Kv33_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dn_Kv33_dt_linearized))*dn_Kv33_dt/dn_Kv33_dt_linearized :
    dt*dn_Kv33_dt) + n_Kv33;
  const double I_Kv33 = g_Kv33*std::pow(n_Kv33, 4.)*(-E_K + v);

  // Expressions for the Kv3.4 component
  const double q10_Kv34 = 3.;
  const double qt_Kv34 = std::pow(q10_Kv34, -3.7 + 0.1*Temp);
  const double mivh_Kv34 = -24.;
  const double mik_Kv34 = 15.4000000000000;
  const double mty0_Kv34 = 0.000128510000000000;
  const double mtvh1_Kv34 = 100.700000000000;
  const double mtk1_Kv34 = 12.9000000000000;
  const double mtvh2_Kv34 = -56.0000000000000;
  const double mtk2_Kv34 = -23.1000000000000;
  const double hiy0_Kv34 = 0.310000000000000;
  const double hiA_Kv34 = 0.690000000000000;
  const double hivh_Kv34 = -5.80200000000000;
  const double hik_Kv34 = 11.2000000000000;
  const double minf_Kv34 = 1.0/(1.0 + std::exp((-11.0 + mivh_Kv34 -
    v)/mik_Kv34));
  const double mtau_Kv34 = 1000.*(11. + v < -35. ? 0.000102675 +
    0.0220402902836158*std::exp(0.0353481795687522*v) : mty0_Kv34 +
    1.0/(std::exp((11.0 + mtvh1_Kv34 + v)/mtk1_Kv34) + std::exp((11.0 +
    mtvh2_Kv34 + v)/mtk2_Kv34)))/qt_Kv34;
  const double hinf_Kv34 = hiy0_Kv34 + hiA_Kv34/(1. + std::exp((11.0 -
    hivh_Kv34 + v)/hik_Kv34));
  const double htau_Kv34 = 1000.*(11. + v > 0. ? 0.0012 +
    0.000487682413465536*std::exp(-0.141*v) : 1.2202e-5 +
    0.012*std::exp(-((1.35685483870968 +
    0.0201612903225806*v)*(1.35685483870968 +
    0.0201612903225806*v))))/qt_Kv34;
  const double dm_Kv34_dt = (-m_Kv34 + minf_Kv34)*allow_change/mtau_Kv34;
  const double dm_Kv34_dt_linearized = -allow_change/mtau_Kv34;
  states[17] = (std::fabs(dm_Kv34_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dm_Kv34_dt_linearized))*dm_Kv34_dt/dm_Kv34_dt_linearized :
    dt*dm_Kv34_dt) + m_Kv34;
  const double dh_Kv34_dt = (-h_Kv34 + hinf_Kv34)*allow_change/htau_Kv34;
  const double dh_Kv34_dt_linearized = -allow_change/htau_Kv34;
  states[18] = (std::fabs(dh_Kv34_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dh_Kv34_dt_linearized))*dh_Kv34_dt/dh_Kv34_dt_linearized :
    dt*dh_Kv34_dt) + h_Kv34;
  const double I_Kv34 = g_Kv34*(m_Kv34*m_Kv34*m_Kv34)*(-E_K + v)*h_Kv34;

  // Expressions for the Kv4.3 component
  const double Q10_Kv43 = std::pow(3., -2.55 + 0.1*Temp);
  const double Aalpha_a_Kv43 = 0.814700000000000;
  const double Kalpha_a_Kv43 = -23.3270800000000;
  const double V0alpha_a_Kv43 = -9.17203000000000;
  const double Abeta_a_Kv43 = 0.165500000000000;
  const double Kbeta_a_Kv43 = 19.4717500000000;
  const double V0beta_a_Kv43 = -18.2791400000000;
  const double Aalpha_b_Kv43 = 0.0368000000000000;
  const double Kalpha_b_Kv43 = 12.8433000000000;
  const double V0alpha_b_Kv43 = -111.332090000000;
  const double Abeta_b_Kv43 = 0.0345000000000000;
  const double Kbeta_b_Kv43 = -8.90123000000000;
  const double V0beta_b_Kv43 = -49.9537000000000;
  const double alpha_a_Kv43 = Aalpha_a_Kv43*Q10_Kv43/(1. +
    std::exp((-V0alpha_a_Kv43 + v)/Kalpha_a_Kv43));
  const double beta_a_Kv43 = Abeta_a_Kv43*Q10_Kv43*std::exp(-(-V0beta_a_Kv43 +
    v)/Kbeta_a_Kv43);
  const double alpha_b_Kv43 = Aalpha_b_Kv43*Q10_Kv43/(1. +
    std::exp((-V0alpha_b_Kv43 + v)/Kalpha_b_Kv43));
  const double beta_b_Kv43 = Abeta_b_Kv43*Q10_Kv43/(1. +
    std::exp((-V0beta_b_Kv43 + v)/Kbeta_b_Kv43));
  const double a_inf_Kv43 = alpha_a_Kv43/(alpha_a_Kv43 + beta_a_Kv43);
  const double tau_a_Kv43 = 1.0/(alpha_a_Kv43 + beta_a_Kv43);
  const double b_inf_Kv43 = alpha_b_Kv43/(alpha_b_Kv43 + beta_b_Kv43);
  const double tau_b_Kv43 = 1.0/(alpha_b_Kv43 + beta_b_Kv43);
  const double da_Kv43_dt = (-a_Kv43 + a_inf_Kv43)*allow_change/tau_a_Kv43;
  const double da_Kv43_dt_linearized = -allow_change/tau_a_Kv43;
  states[19] = (std::fabs(da_Kv43_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*da_Kv43_dt_linearized))*da_Kv43_dt/da_Kv43_dt_linearized :
    dt*da_Kv43_dt) + a_Kv43;
  const double db_Kv43_dt = (-b_Kv43 + b_inf_Kv43)*allow_change/tau_b_Kv43;
  const double db_Kv43_dt_linearized = -allow_change/tau_b_Kv43;
  states[20] = (std::fabs(db_Kv43_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*db_Kv43_dt_linearized))*db_Kv43_dt/db_Kv43_dt_linearized :
    dt*db_Kv43_dt) + b_Kv43;
  const double I_Kv43 = g_Kv43*(a_Kv43*a_Kv43*a_Kv43)*(-E_K + v)*b_Kv43;

  // Expressions for the Kir2.x component
  const double Q10_Kir2x = std::pow(3., -2.0 + 0.1*Temp);
  const double Aalpha_d_Kir2x = 0.132890000000000;
  const double Kalpha_d_Kir2x = -24.3902000000000;
  const double V0alpha_d_Kir2x = -83.9400000000000;
  const double Abeta_d_Kir2x = 0.169940000000000;
  const double Kbeta_d_Kir2x = 35.7140000000000;
  const double V0beta_d_Kir2x = -83.9400000000000;
  const double a_d_Kir2x =
    Aalpha_d_Kir2x*Q10_Kir2x*std::exp((-V0alpha_d_Kir2x + v)/Kalpha_d_Kir2x);
  const double b_d_Kir2x = Abeta_d_Kir2x*Q10_Kir2x*std::exp((-V0beta_d_Kir2x
    + v)/Kbeta_d_Kir2x);
  const double tau_d_Kir2x = 1.0/(a_d_Kir2x + b_d_Kir2x);
  const double d_inf_Kir2x = a_d_Kir2x/(a_d_Kir2x + b_d_Kir2x);
  const double dd_Kir2x_dt = (-d_Kir2x + d_inf_Kir2x)*allow_change/tau_d_Kir2x;
  const double dd_Kir2x_dt_linearized = -allow_change/tau_d_Kir2x;
  states[21] = (std::fabs(dd_Kir2x_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dd_Kir2x_dt_linearized))*dd_Kir2x_dt/dd_Kir2x_dt_linearized :
    dt*dd_Kir2x_dt) + d_Kir2x;
  const double I_Kir2x = g_Kir2x*(-E_K + v)*d_Kir2x;

  // Expressions for the Kca1.1 component
  const double q10_Kca11 = 3.;
  const double qt_Kca11 = std::pow(q10_Kca11, -2.3 + 0.1*Temp);
  const double Qo_Kca11 = 0.730000000000000;
  const double Qc_Kca11 = -0.670000000000000;
  const double k1_Kca11 = 1000.00000000000;
  const double onoffrate_Kca11 = 1.;
  const double Kc_Kca11 = 0.0110000000000000;
  const double Ko_Kca11 = 0.00110000000000000;
  const double pf0_Kca11 = 0.00239000000000000;
  const double pf1_Kca11 = 0.00700000000000000;
  const double pf2_Kca11 = 0.0400000000000000;
  const double pf3_Kca11 = 0.295000000000000;
  const double pf4_Kca11 = 0.557000000000000;
  const double pb0_Kca11 = 3.93600000000000;
  const double pb1_Kca11 = 1.15200000000000;
  const double pb2_Kca11 = 0.659000000000000;
  const double pb3_Kca11 = 0.486000000000000;
  const double pb4_Kca11 = 0.0920000000000000;
  const double c01_Kca11 = 4.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c12_Kca11 = 3.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c23_Kca11 = 2.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c34_Kca11 = cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o01_Kca11 = 4.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o12_Kca11 = 3.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o23_Kca11 = 2.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o34_Kca11 = cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c10_Kca11 = Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c21_Kca11 = 2.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c32_Kca11 = 3.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c43_Kca11 = 4.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o10_Kca11 = Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o21_Kca11 = 2.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o32_Kca11 = 3.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o43_Kca11 = 4.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double alpha_Kca11 = std::exp(FARADAY*Qo_Kca11*v/(R*(273.15 + Temp)));
  const double beta_Kca11 = std::exp(FARADAY*Qc_Kca11*v/(R*(273.15 + Temp)));
  const double f0_Kca11 = pf0_Kca11*alpha_Kca11*qt_Kca11;
  const double f1_Kca11 = pf1_Kca11*alpha_Kca11*qt_Kca11;
  const double f2_Kca11 = pf2_Kca11*alpha_Kca11*qt_Kca11;
  const double f3_Kca11 = pf3_Kca11*alpha_Kca11*qt_Kca11;
  const double f4_Kca11 = pf4_Kca11*alpha_Kca11*qt_Kca11;
  const double b0_Kca11 = pb0_Kca11*beta_Kca11*qt_Kca11;
  const double b1_Kca11 = pb1_Kca11*beta_Kca11*qt_Kca11;
  const double b2_Kca11 = pb2_Kca11*beta_Kca11*qt_Kca11;
  const double b3_Kca11 = pb3_Kca11*beta_Kca11*qt_Kca11;
  const double b4_Kca11 = pb4_Kca11*beta_Kca11*qt_Kca11;
  const double C4_Kca11 = 1. - C0_Kca11 - C1_Kca11 - C2_Kca11 - C3_Kca11 -
    O0_Kca11 - O1_Kca11 - O2_Kca11 - O3_Kca11 - O4_Kca11;
  const double dC0_Kca11_dt = ((-c01_Kca11 - f0_Kca11)*C0_Kca11 +
    C1_Kca11*c10_Kca11 + O0_Kca11*b0_Kca11)*allow_change;
  const double dC0_Kca11_dt_linearized = (-c01_Kca11 - f0_Kca11)*allow_change;
  states[22] = C0_Kca11 + (std::fabs(dC0_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dC0_Kca11_dt_linearized))*dC0_Kca11_dt/dC0_Kca11_dt_linearized
    : dt*dC0_Kca11_dt);
  const double dC1_Kca11_dt = ((-c10_Kca11 - c12_Kca11 - f1_Kca11)*C1_Kca11 +
    C0_Kca11*c01_Kca11 + C2_Kca11*c21_Kca11 + O1_Kca11*b1_Kca11)*allow_change;
  const double dC1_Kca11_dt_linearized = (-c10_Kca11 - c12_Kca11 -
    f1_Kca11)*allow_change;
  states[23] = C1_Kca11 + (std::fabs(dC1_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dC1_Kca11_dt_linearized))*dC1_Kca11_dt/dC1_Kca11_dt_linearized
    : dt*dC1_Kca11_dt);
  const double dC2_Kca11_dt = ((-c21_Kca11 - c23_Kca11 - f2_Kca11)*C2_Kca11 +
    C1_Kca11*c12_Kca11 + C3_Kca11*c32_Kca11 + O2_Kca11*b2_Kca11)*allow_change;
  const double dC2_Kca11_dt_linearized = (-c21_Kca11 - c23_Kca11 -
    f2_Kca11)*allow_change;
  states[24] = C2_Kca11 + (std::fabs(dC2_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dC2_Kca11_dt_linearized))*dC2_Kca11_dt/dC2_Kca11_dt_linearized
    : dt*dC2_Kca11_dt);
  const double dC3_Kca11_dt = ((-c32_Kca11 - c34_Kca11 - f3_Kca11)*C3_Kca11 +
    C2_Kca11*c23_Kca11 + C4_Kca11*c43_Kca11 + O3_Kca11*b3_Kca11)*allow_change;
  const double dC3_Kca11_dt_linearized = (-c32_Kca11 - c34_Kca11 - c43_Kca11 -
    f3_Kca11)*allow_change;
  states[25] = C3_Kca11 + (std::fabs(dC3_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dC3_Kca11_dt_linearized))*dC3_Kca11_dt/dC3_Kca11_dt_linearized
    : dt*dC3_Kca11_dt);
  const double dO0_Kca11_dt = ((-b0_Kca11 - o01_Kca11)*O0_Kca11 +
    C0_Kca11*f0_Kca11 + O1_Kca11*o10_Kca11)*allow_change;
  const double dO0_Kca11_dt_linearized = (-b0_Kca11 - o01_Kca11)*allow_change;
  states[26] = O0_Kca11 + (std::fabs(dO0_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dO0_Kca11_dt_linearized))*dO0_Kca11_dt/dO0_Kca11_dt_linearized
    : dt*dO0_Kca11_dt);
  const double dO1_Kca11_dt = ((-b1_Kca11 - o10_Kca11 - o12_Kca11)*O1_Kca11 +
    C1_Kca11*f1_Kca11 + O0_Kca11*o01_Kca11 + O2_Kca11*o21_Kca11)*allow_change;
  const double dO1_Kca11_dt_linearized = (-b1_Kca11 - o10_Kca11 -
    o12_Kca11)*allow_change;
  states[27] = O1_Kca11 + (std::fabs(dO1_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dO1_Kca11_dt_linearized))*dO1_Kca11_dt/dO1_Kca11_dt_linearized
    : dt*dO1_Kca11_dt);
  const double dO2_Kca11_dt = ((-b2_Kca11 - o21_Kca11 - o23_Kca11)*O2_Kca11 +
    C2_Kca11*f2_Kca11 + O1_Kca11*o12_Kca11 + O3_Kca11*o32_Kca11)*allow_change;
  const double dO2_Kca11_dt_linearized = (-b2_Kca11 - o21_Kca11 -
    o23_Kca11)*allow_change;
  states[28] = O2_Kca11 + (std::fabs(dO2_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dO2_Kca11_dt_linearized))*dO2_Kca11_dt/dO2_Kca11_dt_linearized
    : dt*dO2_Kca11_dt);
  const double dO3_Kca11_dt = ((-b3_Kca11 - o32_Kca11 - o34_Kca11)*O3_Kca11 +
    C3_Kca11*f3_Kca11 + O2_Kca11*o23_Kca11 + O4_Kca11*o43_Kca11)*allow_change;
  const double dO3_Kca11_dt_linearized = (-b3_Kca11 - o32_Kca11 -
    o34_Kca11)*allow_change;
  states[29] = O3_Kca11 + (std::fabs(dO3_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dO3_Kca11_dt_linearized))*dO3_Kca11_dt/dO3_Kca11_dt_linearized
    : dt*dO3_Kca11_dt);
  const double dO4_Kca11_dt = ((-b4_Kca11 - o43_Kca11)*O4_Kca11 +
    C4_Kca11*f4_Kca11 + O3_Kca11*o34_Kca11)*allow_change;
  const double dO4_Kca11_dt_linearized = (-b4_Kca11 - f4_Kca11 -
    o43_Kca11)*allow_change;
  states[30] = O4_Kca11 + (std::fabs(dO4_Kca11_dt_linearized) > 1.0e-8 ?
    (-1.0 +
    std::exp(dt*dO4_Kca11_dt_linearized))*dO4_Kca11_dt/dO4_Kca11_dt_linearized
    : dt*dO4_Kca11_dt);
  const double I_Kca11 = g_Kca11*(-E_K + v)*(O0_Kca11 + O1_Kca11 + O2_Kca11 +
    O3_Kca11 + O4_Kca11);

  // Expressions for the Kca2.2 component
  const double Q10_Kca22 = 3.00000000000000;
  const double tcorr_Kca22 = std::pow(Q10_Kca22, -2.3 + 0.1*Temp);
  const double invc1_Kca22 = 0.0800000000000000;
  const double invc2_Kca22 = 0.0800000000000000;
  const double invc3_Kca22 = 0.200000000000000;
  const double invo1_Kca22 = 1.;
  const double invo2_Kca22 = 0.100000000000000;
  const double diro1_Kca22 = 0.160000000000000;
  const double diro2_Kca22 = 1.20000000000000;
  const double dirc2_Kca22 = 200.000000000000;
  const double dirc3_Kca22 = 160.000000000000;
  const double dirc4_Kca22 = 80.0000000000000;
  const double invc1_t_Kca22 = invc1_Kca22*tcorr_Kca22;
  const double invc2_t_Kca22 = invc2_Kca22*tcorr_Kca22;
  const double invc3_t_Kca22 = invc3_Kca22*tcorr_Kca22;
  const double invo1_t_Kca22 = invo1_Kca22*tcorr_Kca22;
  const double invo2_t_Kca22 = invo2_Kca22*tcorr_Kca22;
  const double diro1_t_Kca22 = diro1_Kca22*tcorr_Kca22;
  const double diro2_t_Kca22 = diro2_Kca22*tcorr_Kca22;
  const double dirc2_t_Kca22 = dirc2_Kca22*tcorr_Kca22;
  const double dirc3_t_Kca22 = dirc3_Kca22*tcorr_Kca22;
  const double dirc4_t_Kca22 = dirc4_Kca22*tcorr_Kca22;
  const double dirc2_t_ca_Kca22 = cai*dirc2_t_Kca22;
  const double dirc3_t_ca_Kca22 = cai*dirc3_t_Kca22;
  const double dirc4_t_ca_Kca22 = cai*dirc4_t_Kca22;
  const double c4_Kca22 = 1. - c1_Kca22 - c2_Kca22 - c3_Kca22 - o1_Kca22 -
    o2_Kca22;
  const double dc1_Kca22_dt = (c2_Kca22*invc1_t_Kca22 -
    c1_Kca22*dirc2_t_ca_Kca22)*allow_change;
  const double dc1_Kca22_dt_linearized = -allow_change*dirc2_t_ca_Kca22;
  states[31] = (std::fabs(dc1_Kca22_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dc1_Kca22_dt_linearized))*dc1_Kca22_dt/dc1_Kca22_dt_linearized
    : dt*dc1_Kca22_dt) + c1_Kca22;
  const double dc2_Kca22_dt = ((-dirc3_t_ca_Kca22 - invc1_t_Kca22)*c2_Kca22 +
    c1_Kca22*dirc2_t_ca_Kca22 + c3_Kca22*invc2_t_Kca22)*allow_change;
  const double dc2_Kca22_dt_linearized = (-dirc3_t_ca_Kca22 -
    invc1_t_Kca22)*allow_change;
  states[32] = (std::fabs(dc2_Kca22_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dc2_Kca22_dt_linearized))*dc2_Kca22_dt/dc2_Kca22_dt_linearized
    : dt*dc2_Kca22_dt) + c2_Kca22;
  const double dc3_Kca22_dt = ((-dirc4_t_ca_Kca22 - diro1_t_Kca22 -
    invc2_t_Kca22)*c3_Kca22 + c2_Kca22*dirc3_t_ca_Kca22 +
    c4_Kca22*invc3_t_Kca22 + invo1_t_Kca22*o1_Kca22)*allow_change;
  const double dc3_Kca22_dt_linearized = (-dirc4_t_ca_Kca22 - diro1_t_Kca22 -
    invc2_t_Kca22 - invc3_t_Kca22)*allow_change;
  states[33] = (std::fabs(dc3_Kca22_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dc3_Kca22_dt_linearized))*dc3_Kca22_dt/dc3_Kca22_dt_linearized
    : dt*dc3_Kca22_dt) + c3_Kca22;
  const double do1_Kca22_dt = (c3_Kca22*diro1_t_Kca22 -
    invo1_t_Kca22*o1_Kca22)*allow_change;
  const double do1_Kca22_dt_linearized = -allow_change*invo1_t_Kca22;
  states[34] = (std::fabs(do1_Kca22_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*do1_Kca22_dt_linearized))*do1_Kca22_dt/do1_Kca22_dt_linearized
    : dt*do1_Kca22_dt) + o1_Kca22;
  const double do2_Kca22_dt = (c4_Kca22*diro2_t_Kca22 -
    invo2_t_Kca22*o2_Kca22)*allow_change;
  const double do2_Kca22_dt_linearized = (-diro2_t_Kca22 -
    invo2_t_Kca22)*allow_change;
  states[35] = (std::fabs(do2_Kca22_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*do2_Kca22_dt_linearized))*do2_Kca22_dt/do2_Kca22_dt_linearized
    : dt*do2_Kca22_dt) + o2_Kca22;
  const double I_Kca22 = g_Kca22*(-E_K + v)*(o1_Kca22 + o2_Kca22);

  // Expressions for the Kca3.1 component
  const double Yvdep_Kca31 = 13.364375107327*std::exp(0.037037037037037*v);
  const double Yconcdep_Kca31 = (cai < 0.01 ? (7.5 - 500.*cai)/(-1. +
    102586.491213197*std::exp(-769.230769230769*cai)) : 0.0545700593112901);
  const double Ybeta_Kca31 = 0.0500000000000000;
  const double Yalpha_Kca31 = Yconcdep_Kca31*Yvdep_Kca31;
  const double dY_Kca31_dt = ((1. - Y_Kca31)*Yalpha_Kca31 -
    Ybeta_Kca31*Y_Kca31)*allow_change;
  const double dY_Kca31_dt_linearized = (-Ybeta_Kca31 -
    Yalpha_Kca31)*allow_change;
  states[36] = (std::fabs(dY_Kca31_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dY_Kca31_dt_linearized))*dY_Kca31_dt/dY_Kca31_dt_linearized :
    dt*dY_Kca31_dt) + Y_Kca31;
  const double I_Kca31 = g_Kca31*(-E_K + v)*Y_Kca31;

  // Expressions for the Cav2.1 component
  const double q10_Cav21 = 3.;
  const double qt_Cav21 = std::pow(q10_Cav21, -2.3 + 0.1*Temp);
  const double vhalfm_Cav21 = -29.4580000000000;
  const double cvm_Cav21 = 8.42900000000000;
  const double vshift_Cav21 = 0.;
  const double minf_Cav21 = 1.0/(1.0 + std::exp((vhalfm_Cav21 + vshift_Cav21 -
    v)/cvm_Cav21));
  const double taum_Cav21 = 1.0*(v >= -40. ? 0.2702 +
    1.1622*std::exp(0.00609050490285645*(26.798 + v)*(-26.798 - v)) :
    0.6923*std::exp(0.00091796007240869*v))/qt_Cav21;
  const double dm_Cav21_dt = (-m_Cav21 + minf_Cav21)*allow_change/taum_Cav21;
  const double dm_Cav21_dt_linearized = -allow_change/taum_Cav21;
  states[37] = (std::fabs(dm_Cav21_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dm_Cav21_dt_linearized))*dm_Cav21_dt/dm_Cav21_dt_linearized :
    dt*dm_Cav21_dt) + m_Cav21;
  const double z_Ca = 2.;
  const double zeta_Ca = FARADAY*z_Ca*v/(R*(273.19 + Temp));
  const double ghk = (std::fabs(1. - std::exp(-zeta_Ca)) < 1.0e-6 ?
    1.0e-6*F*z_Ca*(1. + zeta_Ca/2.)*(cai - cao*std::exp(-zeta_Ca)) :
    1.0e-6*F*z_Ca*(cai - cao*std::exp(-zeta_Ca))*zeta_Ca/(1. -
    std::exp(-zeta_Ca)));
  const double I_Cav21 = 1000.0*g_Cav21*(m_Cav21*m_Cav21*m_Cav21)*ghk;

  // Expressions for the Cav3.1 component
  const double q10_Cav31 = 3.;
  const double qt_Cav31 = std::pow(q10_Cav31, -3.7 + 0.1*Temp);
  const double v0_m_inf_Cav31 = -52.0000000000000;
  const double v0_h_inf_Cav31 = -72.0000000000000;
  const double k_m_inf_Cav31 = -5.00000000000000;
  const double k_h_inf_Cav31 = 7.00000000000000;
  const double C_tau_m_Cav31 = 1.00000000000000;
  const double A_tau_m_Cav31 = 1.00000000000000;
  const double v0_tau_m1_Cav31 = -40.0000000000000;
  const double v0_tau_m2_Cav31 = -102.000000000000;
  const double k_tau_m1_Cav31 = 9.00000000000000;
  const double k_tau_m2_Cav31 = -18.0000000000000;
  const double C_tau_h_Cav31 = 15.0000000000000;
  const double A_tau_h_Cav31 = 1.00000000000000;
  const double v0_tau_h1_Cav31 = -32.0000000000000;
  const double k_tau_h1_Cav31 = 7.00000000000000;
  const double minf_Cav31 = 1.0/(1. + std::exp((-v0_m_inf_Cav31 +
    v)/k_m_inf_Cav31));
  const double hinf_Cav31 = 1.0/(1. + std::exp((-v0_h_inf_Cav31 +
    v)/k_h_inf_Cav31));
  const double taum_Cav31 = (v <= -90. ? 1. : (C_tau_m_Cav31 +
    A_tau_m_Cav31/(std::exp((-v0_tau_m1_Cav31 + v)/k_tau_m1_Cav31) +
    std::exp((-v0_tau_m2_Cav31 + v)/k_tau_m2_Cav31)))/qt_Cav31);
  const double tauh_Cav31 = (C_tau_h_Cav31 +
    A_tau_h_Cav31*std::exp(-(-v0_tau_h1_Cav31 + v)/k_tau_h1_Cav31))/qt_Cav31;
  const double dm_Cav31_dt = (-m_Cav31 + minf_Cav31)*allow_change/taum_Cav31;
  const double dm_Cav31_dt_linearized = -allow_change/taum_Cav31;
  states[38] = (std::fabs(dm_Cav31_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dm_Cav31_dt_linearized))*dm_Cav31_dt/dm_Cav31_dt_linearized :
    dt*dm_Cav31_dt) + m_Cav31;
  const double dh_Cav31_dt = (-h_Cav31 + hinf_Cav31)*allow_change/tauh_Cav31;
  const double dh_Cav31_dt_linearized = -allow_change/tauh_Cav31;
  states[39] = (std::fabs(dh_Cav31_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dh_Cav31_dt_linearized))*dh_Cav31_dt/dh_Cav31_dt_linearized :
    dt*dh_Cav31_dt) + h_Cav31;
  const double I_Cav31 = 1000.0*g_Cav31*(m_Cav31*m_Cav31)*ghk*h_Cav31;

  // Expressions for the Cav3.2 component
  const double phi_m_Cav32 = 6.89864830730607;
  const double phi_h_Cav32 = 3.73719281884655;
  const double m_inf_Cav32 = 1.0/(1. +
    0.000607957605998435*std::exp(-0.135135135135135*v));
  const double h_inf_Cav32 = 1.0/(1. +
    148461.062050627*std::exp(0.139275766016713*v));
  const double tau_m_Cav32 = (1.9 +
    1.0/(22.404093719072*std::exp(0.0840336134453781*v) +
    0.00189854653589182*std::exp(-v/21.)))/phi_m_Cav32;
  const double tau_h_Cav32 = 13.7 + (1942.0 +
    55178666.3027343*std::exp(0.108695652173913*v))/(phi_h_Cav32*(1. +
    30321871936.1927*std::exp(0.27027027027027*v)));
  const double dm_Cav32_dt = (-m_Cav32 + m_inf_Cav32)*allow_change/tau_m_Cav32;
  const double dm_Cav32_dt_linearized = -allow_change/tau_m_Cav32;
  states[40] = (std::fabs(dm_Cav32_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dm_Cav32_dt_linearized))*dm_Cav32_dt/dm_Cav32_dt_linearized :
    dt*dm_Cav32_dt) + m_Cav32;
  const double dh_Cav32_dt = (-h_Cav32 + h_inf_Cav32)*allow_change/tau_h_Cav32;
  const double dh_Cav32_dt_linearized = -allow_change/tau_h_Cav32;
  states[41] = (std::fabs(dh_Cav32_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dh_Cav32_dt_linearized))*dh_Cav32_dt/dh_Cav32_dt_linearized :
    dt*dh_Cav32_dt) + h_Cav32;
  const double E_Ca = R*(273.15 + Temp)*std::log(cao/cai)/(2.*FARADAY);
  const double I_Cav32 = g_Cav32*(m_Cav32*m_Cav32)*(-E_Ca + v)*h_Cav32;

  // Expressions for the Cav3.3 component
  const double q10_Cav33 = 2.30000000000000;
  const double qt_Cav33 = std::pow(q10_Cav33, -2.8 + 0.1*Temp);
  const double gCav3_3bar = 1.00000000000000e-5;
  const double vhalfn_Cav33 = -41.5000000000000;
  const double vhalfl_Cav33 = -69.8000000000000;
  const double kn_Cav33 = 6.20000000000000;
  const double kl_Cav33 = -6.10000000000000;
  const double n_inf_Cav33 = 1.0/(1.0 + std::exp((vhalfn_Cav33 - v)/kn_Cav33));
  const double l_inf_Cav33 = 1.0/(1.0 + std::exp((vhalfl_Cav33 - v)/kl_Cav33));
  const double tau_n_Cav33 = (v > -60. ? 7.2 +
    0.02*std::exp(-0.0680272108843537*v) : 79.5 +
    2.0*std::exp(-0.10752688172043*v))/qt_Cav33;
  const double tau_l_Cav33 = (v > -60. ? 0.875*std::exp(120./41. + v/41.) :
    260.0)/qt_Cav33;
  const double w_Cav33 = FARADAY*z_Ca*v/(R*(273.14 + Temp));
  const double ghk_Cav33 = -FARADAY*z_Ca*(cao -
    cai*std::exp(w_Cav33))*w_Cav33/(-1.0 + std::exp(w_Cav33));
  const double dn_Cav33_dt = (-n_Cav33 + n_inf_Cav33)*allow_change/tau_n_Cav33;
  const double dn_Cav33_dt_linearized = -allow_change/tau_n_Cav33;
  states[42] = (std::fabs(dn_Cav33_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dn_Cav33_dt_linearized))*dn_Cav33_dt/dn_Cav33_dt_linearized :
    dt*dn_Cav33_dt) + n_Cav33;
  const double dl_Cav33_dt = (-l_Cav33 + l_inf_Cav33)*allow_change/tau_l_Cav33;
  const double dl_Cav33_dt_linearized = -allow_change/tau_l_Cav33;
  states[43] = (std::fabs(dl_Cav33_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dl_Cav33_dt_linearized))*dl_Cav33_dt/dl_Cav33_dt_linearized :
    dt*dl_Cav33_dt) + l_Cav33;
  const double I_Cav33 =
    gCav3_3bar*g_Cav33*(n_Cav33*n_Cav33*n_Cav33)*ghk_Cav33*l_Cav33;

  // Expressions for the HCN1 component
  const double q10_HCN1 = 3.00000000000000;
  const double qt_HCN1 = std::pow(q10_HCN1, -3.7 + 0.1*Temp);
  const double ratetau_HCN1 = 1.00000000000000;
  const double ljp_HCN1 = 9.30000000000000;
  const double v_inf_half_noljp_HCN1 = -90.3000000000000;
  const double v_inf_k_HCN1 = 9.67000000000000;
  const double v_tau_const_HCN1 = 0.00180000000000000;
  const double v_tau_half1_noljp_HCN1 = -68.0000000000000;
  const double v_tau_half2_noljp_HCN1 = -68.0000000000000;
  const double v_tau_k1_HCN1 = -22.0000000000000;
  const double v_tau_k2_HCN1 = 7.14000000000000;
  const double v_inf_half_HCN1 = v_inf_half_noljp_HCN1 - ljp_HCN1;
  const double v_tau_half1_HCN1 = v_tau_half1_noljp_HCN1 - ljp_HCN1;
  const double v_tau_half2_HCN1 = v_tau_half2_noljp_HCN1 - ljp_HCN1;
  const double hinf_HCN1 = 1.0/(1. + std::exp((-v_inf_half_HCN1 +
    v)/v_inf_k_HCN1));
  const double tauh_HCN1 =
    ratetau_HCN1/(v_tau_const_HCN1*(std::exp((-v_tau_half1_HCN1 +
    v)/v_tau_k1_HCN1) + std::exp((-v_tau_half2_HCN1 +
    v)/v_tau_k2_HCN1))*qt_HCN1);
  const double dh_HCN1_dt = (-h_HCN1 + hinf_HCN1)*allow_change/tauh_HCN1;
  const double dh_HCN1_dt_linearized = -allow_change/tauh_HCN1;
  states[44] = (std::fabs(dh_HCN1_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dh_HCN1_dt_linearized))*dh_HCN1_dt/dh_HCN1_dt_linearized :
    dt*dh_HCN1_dt) + h_HCN1;
  const double I_HCN1 = g_HCN1*(-E_h + v)*h_HCN1;

  // Expressions for the leak component
  const double I_leak = g_leak*(-E_leak + v);

  // Expressions for the Membrane component
  const double I_tot = I_Cav21 + I_Cav31 + I_Cav32 + I_Cav33 + I_HCN1 +
    I_Kca11 + I_Kca22 + I_Kca31 + I_Kir2x + I_Kv11 + I_Kv15 + I_Kv33 + I_Kv34 +
    I_Kv43 + I_Na + I_leak;
  const double dv_dt = -1.0*I_tot*allow_change/Cm;
  const double dI_Cav32_dv = g_Cav32*(m_Cav32*m_Cav32)*h_Cav32;
  const double dI_Kv11_dv = g_Kv11*std::pow(n_Kv11, 4.);
  const double dI_Kv34_dv = g_Kv34*(m_Kv34*m_Kv34*m_Kv34)*h_Kv34;
  const double dghk_dzeta_Ca = (std::fabs(1. - std::exp(-zeta_Ca)) < 1.0e-6 ?
    5.0e-7*F*z_Ca*(cai - cao*std::exp(-zeta_Ca)) + 1.0e-6*F*cao*z_Ca*(1. +
    zeta_Ca/2.)*std::exp(-zeta_Ca) : 1.0e-6*F*z_Ca*(cai -
    cao*std::exp(-zeta_Ca))/(1. - std::exp(-zeta_Ca)) +
    1.0e-6*F*cao*z_Ca*std::exp(-zeta_Ca)*zeta_Ca/(1. - std::exp(-zeta_Ca)) -
    1.0e-6*F*z_Ca*(cai -
    cao*std::exp(-zeta_Ca))*std::exp(-zeta_Ca)*zeta_Ca/((1. -
    std::exp(-zeta_Ca))*(1. - std::exp(-zeta_Ca))));
  const double dI_Kca11_dv = g_Kca11*(O0_Kca11 + O1_Kca11 + O2_Kca11 +
    O3_Kca11 + O4_Kca11);
  const double dI_Kv43_dv = g_Kv43*(a_Kv43*a_Kv43*a_Kv43)*b_Kv43;
  const double dI_Kca22_dv = g_Kca22*(o1_Kca22 + o2_Kca22);
  const double dghk_Cav33_dw_Cav33 = -FARADAY*z_Ca*(cao -
    cai*std::exp(w_Cav33))/(-1.0 + std::exp(w_Cav33)) +
    FARADAY*cai*z_Ca*std::exp(w_Cav33)*w_Cav33/(-1.0 + std::exp(w_Cav33)) +
    FARADAY*z_Ca*(cao -
    cai*std::exp(w_Cav33))*std::exp(w_Cav33)*w_Cav33/((-1.0 +
    std::exp(w_Cav33))*(-1.0 + std::exp(w_Cav33)));
  const double dI_Cav33_dghk_Cav33 =
    gCav3_3bar*g_Cav33*(n_Cav33*n_Cav33*n_Cav33)*l_Cav33;
  const double dI_Cav21_dghk = 1000.0*g_Cav21*(m_Cav21*m_Cav21*m_Cav21);
  const double dI_Cav31_dghk = 1000.0*g_Cav31*(m_Cav31*m_Cav31)*h_Cav31;
  const double dI_Kv15_dv = g_Kv15*(m_Kv15*m_Kv15*m_Kv15)*(0.1 + 1.0/(1.0 +
    3.17036319488807*std::exp(-0.0769230769230769*v)))*n_Kv15*u_Kv15 +
    0.243874091914467*g_Kv15*(m_Kv15*m_Kv15*m_Kv15)*(-E_K +
    v)*std::exp(-0.0769230769230769*v)*n_Kv15*u_Kv15/((1.0 +
    3.17036319488807*std::exp(-0.0769230769230769*v))*(1.0 +
    3.17036319488807*std::exp(-0.0769230769230769*v)));
  const double dzeta_Ca_dv = FARADAY*z_Ca/(R*(273.19 + Temp));
  const double dw_Cav33_dv = FARADAY*z_Ca/(R*(273.14 + Temp));
  const double dI_Kv33_dv = g_Kv33*std::pow(n_Kv33, 4.);
  const double dv_dt_linearized = -1.0*(g_leak + g_HCN1*h_HCN1 +
    g_Kca31*Y_Kca31 + g_Kir2x*d_Kir2x + g_Na*O_Na +
    dI_Cav21_dghk*dghk_dzeta_Ca*dzeta_Ca_dv +
    dI_Cav31_dghk*dghk_dzeta_Ca*dzeta_Ca_dv +
    dI_Cav33_dghk_Cav33*dghk_Cav33_dw_Cav33*dw_Cav33_dv + dI_Cav32_dv +
    dI_Kca11_dv + dI_Kca22_dv + dI_Kv11_dv + dI_Kv15_dv + dI_Kv33_dv +
    dI_Kv34_dv + dI_Kv43_dv)*allow_change/Cm;
  states[45] = (std::fabs(dv_dt_linearized) > 1.0e-8 ? (-1.0 +
    std::exp(dt*dv_dt_linearized))*dv_dt/dv_dt_linearized : dt*dv_dt) + v;
}

// Compute a forward step using the explicit Euler algorithm to the
// Masoli2015 ODE
void forward_explicit_euler(double* states, const double t, const double dt,
  const double* parameters)
{

  // Assign states
  const double C1_Na = states[0];
  const double C2_Na = states[1];
  const double C3_Na = states[2];
  const double C4_Na = states[3];
  const double C5_Na = states[4];
  const double O_Na = states[5];
  const double B_Na = states[6];
  const double I1_Na = states[7];
  const double I2_Na = states[8];
  const double I3_Na = states[9];
  const double I4_Na = states[10];
  const double I5_Na = states[11];
  const double n_Kv11 = states[12];
  const double m_Kv15 = states[13];
  const double n_Kv15 = states[14];
  const double u_Kv15 = states[15];
  const double n_Kv33 = states[16];
  const double m_Kv34 = states[17];
  const double h_Kv34 = states[18];
  const double a_Kv43 = states[19];
  const double b_Kv43 = states[20];
  const double d_Kir2x = states[21];
  const double C0_Kca11 = states[22];
  const double C1_Kca11 = states[23];
  const double C2_Kca11 = states[24];
  const double C3_Kca11 = states[25];
  const double O0_Kca11 = states[26];
  const double O1_Kca11 = states[27];
  const double O2_Kca11 = states[28];
  const double O3_Kca11 = states[29];
  const double O4_Kca11 = states[30];
  const double c1_Kca22 = states[31];
  const double c2_Kca22 = states[32];
  const double c3_Kca22 = states[33];
  const double o1_Kca22 = states[34];
  const double o2_Kca22 = states[35];
  const double Y_Kca31 = states[36];
  const double m_Cav21 = states[37];
  const double m_Cav31 = states[38];
  const double h_Cav31 = states[39];
  const double m_Cav32 = states[40];
  const double h_Cav32 = states[41];
  const double n_Cav33 = states[42];
  const double l_Cav33 = states[43];
  const double h_HCN1 = states[44];
  const double v = states[45];

  // Assign parameters
  const double Cm = parameters[0];
  const double F = parameters[1];
  const double FARADAY = parameters[2];
  const double R = parameters[3];
  const double Tdiff = parameters[4];
  const double Temp = parameters[5];
  const double cai = parameters[6];
  const double cao = parameters[7];
  const double E_Na = parameters[8];
  const double g_Na = parameters[9];
  const double E_K = parameters[10];
  const double g_Kv11 = parameters[11];
  const double g_Kv15 = parameters[12];
  const double g_Kv33 = parameters[13];
  const double g_Kv34 = parameters[14];
  const double g_Kv43 = parameters[15];
  const double g_Kir2x = parameters[16];
  const double g_Kca11 = parameters[17];
  const double g_Kca22 = parameters[18];
  const double g_Kca31 = parameters[19];
  const double g_Cav21 = parameters[20];
  const double g_Cav31 = parameters[21];
  const double g_Cav32 = parameters[22];
  const double g_Cav33 = parameters[23];
  const double E_h = parameters[24];
  const double g_HCN1 = parameters[25];
  const double E_leak = parameters[26];
  const double g_leak = parameters[27];

  // Expressions for the Nav1.6 component
  const double allow_change = (t >= Tdiff ? 1. : 0.);
  const double q10_Na = 3.;
  const double qt_Na = std::pow(q10_Na, -2.2 + 0.1*Temp);
  const double Con_Na = 0.00500000000000000;
  const double Coff_Na = 0.500000000000000;
  const double Oon_Na = 0.750000000000000;
  const double Ooff_Na = 0.00500000000000000;
  const double alpha_Na = 150.;
  const double beta_Na = 3.;
  const double gamma_Na = 150.;
  const double delta_Na = 40.;
  const double epsilon_Na = 1.75000000000000;
  const double zeta_Na = 0.0300000000000000;
  const double x1_Na = 20.;
  const double x2_Na = -20.;
  const double x3_Na = 1000000000000.00;
  const double x4_Na = -1000000000000.00;
  const double x5_Na = 1000000000000.00;
  const double x6_Na = -25.;
  const double alfac_Na = std::pow(Oon_Na/Con_Na, 0.25);
  const double btfac_Na = std::pow(Ooff_Na/Coff_Na, 0.25);
  const double f01_Na = 4.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f02_Na = 3.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f03_Na = 2.*alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f04_Na = alpha_Na*std::exp(v/x1_Na)*qt_Na;
  const double f0O_Na = gamma_Na*std::exp(v/x3_Na)*qt_Na;
  const double fip_Na = epsilon_Na*std::exp(v/x5_Na)*qt_Na;
  const double f11_Na = 4.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f12_Na = 3.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f13_Na = 2.*alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f14_Na = alpha_Na*alfac_Na*std::exp(v/x1_Na)*qt_Na;
  const double f1n_Na = gamma_Na*std::exp(v/x3_Na)*qt_Na;
  const double fi1_Na = Con_Na*qt_Na;
  const double fi2_Na = Con_Na*alfac_Na*qt_Na;
  const double fi3_Na = Con_Na*(alfac_Na*alfac_Na)*qt_Na;
  const double fi4_Na = Con_Na*(alfac_Na*alfac_Na*alfac_Na)*qt_Na;
  const double fi5_Na = Con_Na*std::pow(alfac_Na, 4.)*qt_Na;
  const double fin_Na = Oon_Na*qt_Na;
  const double b01_Na = beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b02_Na = 2.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b03_Na = 3.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b04_Na = 4.*beta_Na*std::exp(v/x2_Na)*qt_Na;
  const double b0O_Na = delta_Na*std::exp(v/x4_Na)*qt_Na;
  const double bip_Na = zeta_Na*std::exp(v/x6_Na)*qt_Na;
  const double b11_Na = beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b12_Na = 2.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b13_Na = 3.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b14_Na = 4.*beta_Na*btfac_Na*std::exp(v/x2_Na)*qt_Na;
  const double b1n_Na = delta_Na*std::exp(v/x4_Na)*qt_Na;
  const double bi1_Na = Coff_Na*qt_Na;
  const double bi2_Na = Coff_Na*btfac_Na*qt_Na;
  const double bi3_Na = Coff_Na*(btfac_Na*btfac_Na)*qt_Na;
  const double bi4_Na = Coff_Na*(btfac_Na*btfac_Na*btfac_Na)*qt_Na;
  const double bi5_Na = Coff_Na*std::pow(btfac_Na, 4.)*qt_Na;
  const double bin_Na = Ooff_Na*qt_Na;
  const double I6_Na = 1. - B_Na - C1_Na - C2_Na - C3_Na - C4_Na - C5_Na -
    I1_Na - I2_Na - I3_Na - I4_Na - I5_Na - O_Na;
  const double dC1_Na_dt = ((-f01_Na - fi1_Na)*C1_Na + C2_Na*b01_Na +
    I1_Na*bi1_Na)*allow_change;
  states[0] = dt*dC1_Na_dt + C1_Na;
  const double dC2_Na_dt = ((-b01_Na - f02_Na - fi2_Na)*C2_Na + C1_Na*f01_Na
    + C3_Na*b02_Na + I2_Na*bi2_Na)*allow_change;
  states[1] = dt*dC2_Na_dt + C2_Na;
  const double dC3_Na_dt = ((-b02_Na - f03_Na - fi3_Na)*C3_Na + C2_Na*f02_Na
    + C4_Na*b03_Na + I3_Na*bi3_Na)*allow_change;
  states[2] = dt*dC3_Na_dt + C3_Na;
  const double dC4_Na_dt = ((-b03_Na - f04_Na - fi4_Na)*C4_Na + C3_Na*f03_Na
    + C5_Na*b04_Na + I4_Na*bi4_Na)*allow_change;
  states[3] = dt*dC4_Na_dt + C4_Na;
  const double dC5_Na_dt = ((-b04_Na - f0O_Na - fi5_Na)*C5_Na + C4_Na*f04_Na
    + I5_Na*bi5_Na + O_Na*b0O_Na)*allow_change;
  states[4] = dt*dC5_Na_dt + C5_Na;
  const double dO_Na_dt = ((-b0O_Na - fin_Na - fip_Na)*O_Na + B_Na*bip_Na +
    C5_Na*f0O_Na + I6_Na*bin_Na)*allow_change;
  states[5] = dt*dO_Na_dt + O_Na;
  const double dB_Na_dt = (O_Na*fip_Na - B_Na*bip_Na)*allow_change;
  states[6] = dt*dB_Na_dt + B_Na;
  const double dI1_Na_dt = ((-bi1_Na - f11_Na)*I1_Na + C1_Na*fi1_Na +
    I2_Na*b11_Na)*allow_change;
  states[7] = dt*dI1_Na_dt + I1_Na;
  const double dI2_Na_dt = ((-b11_Na - bi2_Na - f12_Na)*I2_Na + C2_Na*fi2_Na
    + I1_Na*f11_Na + I3_Na*b12_Na)*allow_change;
  states[8] = dt*dI2_Na_dt + I2_Na;
  const double dI3_Na_dt = ((-b12_Na - bi3_Na - f13_Na)*I3_Na + C3_Na*fi3_Na
    + I2_Na*f12_Na + I4_Na*b13_Na)*allow_change;
  states[9] = dt*dI3_Na_dt + I3_Na;
  const double dI4_Na_dt = ((-b13_Na - bi4_Na - f14_Na)*I4_Na + C4_Na*fi4_Na
    + I3_Na*f13_Na + I5_Na*b14_Na)*allow_change;
  states[10] = dt*dI4_Na_dt + I4_Na;
  const double dI5_Na_dt = ((-b14_Na - bi5_Na - f1n_Na)*I5_Na + C5_Na*fi5_Na
    + I4_Na*f14_Na + I6_Na*b1n_Na)*allow_change;
  states[11] = dt*dI5_Na_dt + I5_Na;
  const double I_Na = g_Na*(-E_Na + v)*O_Na;

  // Expressions for the Kv1.1 component
  const double q10_Kv11 = 2.70000000000000;
  const double ca_Kv11 = 0.128890000000000;
  const double cva_Kv11 = 45.;
  const double cka_Kv11 = -33.9087700000000;
  const double cb_Kv11 = 0.128890000000000;
  const double cvb_Kv11 = 45.;
  const double ckb_Kv11 = 12.4210100000000;
  const double qt_Kv11 = std::pow(q10_Kv11, -2.2 + 0.1*Temp);
  const double alphan_Kv11 = ca_Kv11*std::exp((-cva_Kv11 - v)/cka_Kv11);
  const double betan_Kv11 = cb_Kv11*std::exp((-cvb_Kv11 - v)/ckb_Kv11);
  const double ninf_Kv11 = alphan_Kv11/(alphan_Kv11 + betan_Kv11);
  const double taun_Kv11 = 1./((alphan_Kv11 + betan_Kv11)*qt_Kv11);
  const double dn_Kv11_dt = (-n_Kv11 + ninf_Kv11)*allow_change/taun_Kv11;
  states[12] = dt*dn_Kv11_dt + n_Kv11;
  const double I_Kv11 = g_Kv11*std::pow(n_Kv11, 4.)*(-E_K + v);

  // Expressions for the Kv1.5 component
  const double q10_Kv15 = std::pow(2.2, -3.7 + 0.1*Temp);
  const double am_Kv15 =
    0.65*q10_Kv15/(1.66275285675459*std::exp(-0.0169491525423729*v) +
    0.308365167896581*std::exp(-0.117647058823529*v));
  const double bm_Kv15 = 0.65*q10_Kv15/(2.5 +
    124.403387639703*std::exp(0.0588235294117647*v));
  const double mtau_Kv15 = 1./(3.*(am_Kv15 + bm_Kv15));
  const double minf_Kv15 = 1.0/(1.0 +
    0.0425851362887876*std::exp(-0.104166666666667*v));
  const double an_Kv15 = 0.001*q10_Kv15/(2.4 +
    3.43809189356453*std::exp(-0.0128205128205128*v));
  const double bn_Kv15 = 2.75364493497472e-8*std::exp(0.0625*v)*q10_Kv15;
  const double ntau_Kv15 = 0.333333333333333/(an_Kv15 + bn_Kv15);
  const double ninf_Kv15 = 0.25 + 1.0/(1.35 +
    1.64872127070013*std::exp(0.0714285714285714*v));
  const double uinf_Kv15 = 0.1 + 1.0/(1.1 +
    1.64872127070013*std::exp(0.0714285714285714*v));
  const double utau_Kv15 = 6800.00000000000;
  const double dm_Kv15_dt = (-m_Kv15 + minf_Kv15)*allow_change/mtau_Kv15;
  states[13] = dt*dm_Kv15_dt + m_Kv15;
  const double dn_Kv15_dt = (-n_Kv15 + ninf_Kv15)*allow_change/ntau_Kv15;
  states[14] = dt*dn_Kv15_dt + n_Kv15;
  const double du_Kv15_dt = (-u_Kv15 + uinf_Kv15)*allow_change/utau_Kv15;
  states[15] = dt*du_Kv15_dt + u_Kv15;
  const double I_Kv15 = g_Kv15*(m_Kv15*m_Kv15*m_Kv15)*(0.1 + 1.0/(1.0 +
    3.17036319488807*std::exp(-0.0769230769230769*v)))*(-E_K +
    v)*n_Kv15*u_Kv15;

  // Expressions for the Kv3.3 component
  const double q10_Kv33 = 2.70000000000000;
  const double qt_Kv33 = std::pow(q10_Kv33, -2.2 + 0.1*Temp);
  const double ca_Kv33 = 0.220000000000000;
  const double cva_Kv33 = 16.;
  const double cka_Kv33 = -26.5000000000000;
  const double cb_Kv33 = 0.220000000000000;
  const double cvb_Kv33 = 16.;
  const double ckb_Kv33 = 26.5000000000000;
  const double alpha_Kv33 = ca_Kv33*std::exp((-cva_Kv33 - v)/cka_Kv33)*qt_Kv33;
  const double beta_Kv33 = cb_Kv33*std::exp((-cvb_Kv33 - v)/ckb_Kv33)*qt_Kv33;
  const double dn_Kv33_dt = ((1. - n_Kv33)*alpha_Kv33 -
    beta_Kv33*n_Kv33)*allow_change;
  states[16] = dt*dn_Kv33_dt + n_Kv33;
  const double I_Kv33 = g_Kv33*std::pow(n_Kv33, 4.)*(-E_K + v);

  // Expressions for the Kv3.4 component
  const double q10_Kv34 = 3.;
  const double qt_Kv34 = std::pow(q10_Kv34, -3.7 + 0.1*Temp);
  const double mivh_Kv34 = -24.;
  const double mik_Kv34 = 15.4000000000000;
  const double mty0_Kv34 = 0.000128510000000000;
  const double mtvh1_Kv34 = 100.700000000000;
  const double mtk1_Kv34 = 12.9000000000000;
  const double mtvh2_Kv34 = -56.0000000000000;
  const double mtk2_Kv34 = -23.1000000000000;
  const double hiy0_Kv34 = 0.310000000000000;
  const double hiA_Kv34 = 0.690000000000000;
  const double hivh_Kv34 = -5.80200000000000;
  const double hik_Kv34 = 11.2000000000000;
  const double minf_Kv34 = 1.0/(1.0 + std::exp((-11.0 + mivh_Kv34 -
    v)/mik_Kv34));
  const double mtau_Kv34 = 1000.*(11. + v < -35. ? 0.000102675 +
    0.0220402902836158*std::exp(0.0353481795687522*v) : mty0_Kv34 +
    1.0/(std::exp((11.0 + mtvh1_Kv34 + v)/mtk1_Kv34) + std::exp((11.0 +
    mtvh2_Kv34 + v)/mtk2_Kv34)))/qt_Kv34;
  const double hinf_Kv34 = hiy0_Kv34 + hiA_Kv34/(1. + std::exp((11.0 -
    hivh_Kv34 + v)/hik_Kv34));
  const double htau_Kv34 = 1000.*(11. + v > 0. ? 0.0012 +
    0.000487682413465536*std::exp(-0.141*v) : 1.2202e-5 +
    0.012*std::exp(-((1.35685483870968 +
    0.0201612903225806*v)*(1.35685483870968 +
    0.0201612903225806*v))))/qt_Kv34;
  const double dm_Kv34_dt = (-m_Kv34 + minf_Kv34)*allow_change/mtau_Kv34;
  states[17] = dt*dm_Kv34_dt + m_Kv34;
  const double dh_Kv34_dt = (-h_Kv34 + hinf_Kv34)*allow_change/htau_Kv34;
  states[18] = dt*dh_Kv34_dt + h_Kv34;
  const double I_Kv34 = g_Kv34*(m_Kv34*m_Kv34*m_Kv34)*(-E_K + v)*h_Kv34;

  // Expressions for the Kv4.3 component
  const double Q10_Kv43 = std::pow(3., -2.55 + 0.1*Temp);
  const double Aalpha_a_Kv43 = 0.814700000000000;
  const double Kalpha_a_Kv43 = -23.3270800000000;
  const double V0alpha_a_Kv43 = -9.17203000000000;
  const double Abeta_a_Kv43 = 0.165500000000000;
  const double Kbeta_a_Kv43 = 19.4717500000000;
  const double V0beta_a_Kv43 = -18.2791400000000;
  const double Aalpha_b_Kv43 = 0.0368000000000000;
  const double Kalpha_b_Kv43 = 12.8433000000000;
  const double V0alpha_b_Kv43 = -111.332090000000;
  const double Abeta_b_Kv43 = 0.0345000000000000;
  const double Kbeta_b_Kv43 = -8.90123000000000;
  const double V0beta_b_Kv43 = -49.9537000000000;
  const double alpha_a_Kv43 = Aalpha_a_Kv43*Q10_Kv43/(1. +
    std::exp((-V0alpha_a_Kv43 + v)/Kalpha_a_Kv43));
  const double beta_a_Kv43 = Abeta_a_Kv43*Q10_Kv43*std::exp(-(-V0beta_a_Kv43 +
    v)/Kbeta_a_Kv43);
  const double alpha_b_Kv43 = Aalpha_b_Kv43*Q10_Kv43/(1. +
    std::exp((-V0alpha_b_Kv43 + v)/Kalpha_b_Kv43));
  const double beta_b_Kv43 = Abeta_b_Kv43*Q10_Kv43/(1. +
    std::exp((-V0beta_b_Kv43 + v)/Kbeta_b_Kv43));
  const double a_inf_Kv43 = alpha_a_Kv43/(alpha_a_Kv43 + beta_a_Kv43);
  const double tau_a_Kv43 = 1.0/(alpha_a_Kv43 + beta_a_Kv43);
  const double b_inf_Kv43 = alpha_b_Kv43/(alpha_b_Kv43 + beta_b_Kv43);
  const double tau_b_Kv43 = 1.0/(alpha_b_Kv43 + beta_b_Kv43);
  const double da_Kv43_dt = (-a_Kv43 + a_inf_Kv43)*allow_change/tau_a_Kv43;
  states[19] = dt*da_Kv43_dt + a_Kv43;
  const double db_Kv43_dt = (-b_Kv43 + b_inf_Kv43)*allow_change/tau_b_Kv43;
  states[20] = dt*db_Kv43_dt + b_Kv43;
  const double I_Kv43 = g_Kv43*(a_Kv43*a_Kv43*a_Kv43)*(-E_K + v)*b_Kv43;

  // Expressions for the Kir2.x component
  const double Q10_Kir2x = std::pow(3., -2.0 + 0.1*Temp);
  const double Aalpha_d_Kir2x = 0.132890000000000;
  const double Kalpha_d_Kir2x = -24.3902000000000;
  const double V0alpha_d_Kir2x = -83.9400000000000;
  const double Abeta_d_Kir2x = 0.169940000000000;
  const double Kbeta_d_Kir2x = 35.7140000000000;
  const double V0beta_d_Kir2x = -83.9400000000000;
  const double a_d_Kir2x =
    Aalpha_d_Kir2x*Q10_Kir2x*std::exp((-V0alpha_d_Kir2x + v)/Kalpha_d_Kir2x);
  const double b_d_Kir2x = Abeta_d_Kir2x*Q10_Kir2x*std::exp((-V0beta_d_Kir2x
    + v)/Kbeta_d_Kir2x);
  const double tau_d_Kir2x = 1.0/(a_d_Kir2x + b_d_Kir2x);
  const double d_inf_Kir2x = a_d_Kir2x/(a_d_Kir2x + b_d_Kir2x);
  const double dd_Kir2x_dt = (-d_Kir2x + d_inf_Kir2x)*allow_change/tau_d_Kir2x;
  states[21] = dt*dd_Kir2x_dt + d_Kir2x;
  const double I_Kir2x = g_Kir2x*(-E_K + v)*d_Kir2x;

  // Expressions for the Kca1.1 component
  const double q10_Kca11 = 3.;
  const double qt_Kca11 = std::pow(q10_Kca11, -2.3 + 0.1*Temp);
  const double Qo_Kca11 = 0.730000000000000;
  const double Qc_Kca11 = -0.670000000000000;
  const double k1_Kca11 = 1000.00000000000;
  const double onoffrate_Kca11 = 1.;
  const double Kc_Kca11 = 0.0110000000000000;
  const double Ko_Kca11 = 0.00110000000000000;
  const double pf0_Kca11 = 0.00239000000000000;
  const double pf1_Kca11 = 0.00700000000000000;
  const double pf2_Kca11 = 0.0400000000000000;
  const double pf3_Kca11 = 0.295000000000000;
  const double pf4_Kca11 = 0.557000000000000;
  const double pb0_Kca11 = 3.93600000000000;
  const double pb1_Kca11 = 1.15200000000000;
  const double pb2_Kca11 = 0.659000000000000;
  const double pb3_Kca11 = 0.486000000000000;
  const double pb4_Kca11 = 0.0920000000000000;
  const double c01_Kca11 = 4.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c12_Kca11 = 3.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c23_Kca11 = 2.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c34_Kca11 = cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o01_Kca11 = 4.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o12_Kca11 = 3.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o23_Kca11 = 2.*cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o34_Kca11 = cai*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c10_Kca11 = Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c21_Kca11 = 2.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c32_Kca11 = 3.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double c43_Kca11 = 4.*Kc_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o10_Kca11 = Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o21_Kca11 = 2.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o32_Kca11 = 3.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double o43_Kca11 = 4.*Ko_Kca11*k1_Kca11*onoffrate_Kca11*qt_Kca11;
  const double alpha_Kca11 = std::exp(FARADAY*Qo_Kca11*v/(R*(273.15 + Temp)));
  const double beta_Kca11 = std::exp(FARADAY*Qc_Kca11*v/(R*(273.15 + Temp)));
  const double f0_Kca11 = pf0_Kca11*alpha_Kca11*qt_Kca11;
  const double f1_Kca11 = pf1_Kca11*alpha_Kca11*qt_Kca11;
  const double f2_Kca11 = pf2_Kca11*alpha_Kca11*qt_Kca11;
  const double f3_Kca11 = pf3_Kca11*alpha_Kca11*qt_Kca11;
  const double f4_Kca11 = pf4_Kca11*alpha_Kca11*qt_Kca11;
  const double b0_Kca11 = pb0_Kca11*beta_Kca11*qt_Kca11;
  const double b1_Kca11 = pb1_Kca11*beta_Kca11*qt_Kca11;
  const double b2_Kca11 = pb2_Kca11*beta_Kca11*qt_Kca11;
  const double b3_Kca11 = pb3_Kca11*beta_Kca11*qt_Kca11;
  const double b4_Kca11 = pb4_Kca11*beta_Kca11*qt_Kca11;
  const double C4_Kca11 = 1. - C0_Kca11 - C1_Kca11 - C2_Kca11 - C3_Kca11 -
    O0_Kca11 - O1_Kca11 - O2_Kca11 - O3_Kca11 - O4_Kca11;
  const double dC0_Kca11_dt = ((-c01_Kca11 - f0_Kca11)*C0_Kca11 +
    C1_Kca11*c10_Kca11 + O0_Kca11*b0_Kca11)*allow_change;
  states[22] = dt*dC0_Kca11_dt + C0_Kca11;
  const double dC1_Kca11_dt = ((-c10_Kca11 - c12_Kca11 - f1_Kca11)*C1_Kca11 +
    C0_Kca11*c01_Kca11 + C2_Kca11*c21_Kca11 + O1_Kca11*b1_Kca11)*allow_change;
  states[23] = dt*dC1_Kca11_dt + C1_Kca11;
  const double dC2_Kca11_dt = ((-c21_Kca11 - c23_Kca11 - f2_Kca11)*C2_Kca11 +
    C1_Kca11*c12_Kca11 + C3_Kca11*c32_Kca11 + O2_Kca11*b2_Kca11)*allow_change;
  states[24] = dt*dC2_Kca11_dt + C2_Kca11;
  const double dC3_Kca11_dt = ((-c32_Kca11 - c34_Kca11 - f3_Kca11)*C3_Kca11 +
    C2_Kca11*c23_Kca11 + C4_Kca11*c43_Kca11 + O3_Kca11*b3_Kca11)*allow_change;
  states[25] = dt*dC3_Kca11_dt + C3_Kca11;
  const double dO0_Kca11_dt = ((-b0_Kca11 - o01_Kca11)*O0_Kca11 +
    C0_Kca11*f0_Kca11 + O1_Kca11*o10_Kca11)*allow_change;
  states[26] = dt*dO0_Kca11_dt + O0_Kca11;
  const double dO1_Kca11_dt = ((-b1_Kca11 - o10_Kca11 - o12_Kca11)*O1_Kca11 +
    C1_Kca11*f1_Kca11 + O0_Kca11*o01_Kca11 + O2_Kca11*o21_Kca11)*allow_change;
  states[27] = dt*dO1_Kca11_dt + O1_Kca11;
  const double dO2_Kca11_dt = ((-b2_Kca11 - o21_Kca11 - o23_Kca11)*O2_Kca11 +
    C2_Kca11*f2_Kca11 + O1_Kca11*o12_Kca11 + O3_Kca11*o32_Kca11)*allow_change;
  states[28] = dt*dO2_Kca11_dt + O2_Kca11;
  const double dO3_Kca11_dt = ((-b3_Kca11 - o32_Kca11 - o34_Kca11)*O3_Kca11 +
    C3_Kca11*f3_Kca11 + O2_Kca11*o23_Kca11 + O4_Kca11*o43_Kca11)*allow_change;
  states[29] = dt*dO3_Kca11_dt + O3_Kca11;
  const double dO4_Kca11_dt = ((-b4_Kca11 - o43_Kca11)*O4_Kca11 +
    C4_Kca11*f4_Kca11 + O3_Kca11*o34_Kca11)*allow_change;
  states[30] = dt*dO4_Kca11_dt + O4_Kca11;
  const double I_Kca11 = g_Kca11*(-E_K + v)*(O0_Kca11 + O1_Kca11 + O2_Kca11 +
    O3_Kca11 + O4_Kca11);

  // Expressions for the Kca2.2 component
  const double Q10_Kca22 = 3.00000000000000;
  const double tcorr_Kca22 = std::pow(Q10_Kca22, -2.3 + 0.1*Temp);
  const double invc1_Kca22 = 0.0800000000000000;
  const double invc2_Kca22 = 0.0800000000000000;
  const double invc3_Kca22 = 0.200000000000000;
  const double invo1_Kca22 = 1.;
  const double invo2_Kca22 = 0.100000000000000;
  const double diro1_Kca22 = 0.160000000000000;
  const double diro2_Kca22 = 1.20000000000000;
  const double dirc2_Kca22 = 200.000000000000;
  const double dirc3_Kca22 = 160.000000000000;
  const double dirc4_Kca22 = 80.0000000000000;
  const double invc1_t_Kca22 = invc1_Kca22*tcorr_Kca22;
  const double invc2_t_Kca22 = invc2_Kca22*tcorr_Kca22;
  const double invc3_t_Kca22 = invc3_Kca22*tcorr_Kca22;
  const double invo1_t_Kca22 = invo1_Kca22*tcorr_Kca22;
  const double invo2_t_Kca22 = invo2_Kca22*tcorr_Kca22;
  const double diro1_t_Kca22 = diro1_Kca22*tcorr_Kca22;
  const double diro2_t_Kca22 = diro2_Kca22*tcorr_Kca22;
  const double dirc2_t_Kca22 = dirc2_Kca22*tcorr_Kca22;
  const double dirc3_t_Kca22 = dirc3_Kca22*tcorr_Kca22;
  const double dirc4_t_Kca22 = dirc4_Kca22*tcorr_Kca22;
  const double dirc2_t_ca_Kca22 = cai*dirc2_t_Kca22;
  const double dirc3_t_ca_Kca22 = cai*dirc3_t_Kca22;
  const double dirc4_t_ca_Kca22 = cai*dirc4_t_Kca22;
  const double c4_Kca22 = 1. - c1_Kca22 - c2_Kca22 - c3_Kca22 - o1_Kca22 -
    o2_Kca22;
  const double dc1_Kca22_dt = (c2_Kca22*invc1_t_Kca22 -
    c1_Kca22*dirc2_t_ca_Kca22)*allow_change;
  states[31] = dt*dc1_Kca22_dt + c1_Kca22;
  const double dc2_Kca22_dt = ((-dirc3_t_ca_Kca22 - invc1_t_Kca22)*c2_Kca22 +
    c1_Kca22*dirc2_t_ca_Kca22 + c3_Kca22*invc2_t_Kca22)*allow_change;
  states[32] = dt*dc2_Kca22_dt + c2_Kca22;
  const double dc3_Kca22_dt = ((-dirc4_t_ca_Kca22 - diro1_t_Kca22 -
    invc2_t_Kca22)*c3_Kca22 + c2_Kca22*dirc3_t_ca_Kca22 +
    c4_Kca22*invc3_t_Kca22 + invo1_t_Kca22*o1_Kca22)*allow_change;
  states[33] = dt*dc3_Kca22_dt + c3_Kca22;
  const double do1_Kca22_dt = (c3_Kca22*diro1_t_Kca22 -
    invo1_t_Kca22*o1_Kca22)*allow_change;
  states[34] = dt*do1_Kca22_dt + o1_Kca22;
  const double do2_Kca22_dt = (c4_Kca22*diro2_t_Kca22 -
    invo2_t_Kca22*o2_Kca22)*allow_change;
  states[35] = dt*do2_Kca22_dt + o2_Kca22;
  const double I_Kca22 = g_Kca22*(-E_K + v)*(o1_Kca22 + o2_Kca22);

  // Expressions for the Kca3.1 component
  const double Yvdep_Kca31 = 13.364375107327*std::exp(0.037037037037037*v);
  const double Yconcdep_Kca31 = (cai < 0.01 ? (7.5 - 500.*cai)/(-1. +
    102586.491213197*std::exp(-769.230769230769*cai)) : 0.0545700593112901);
  const double Ybeta_Kca31 = 0.0500000000000000;
  const double Yalpha_Kca31 = Yconcdep_Kca31*Yvdep_Kca31;
  const double dY_Kca31_dt = ((1. - Y_Kca31)*Yalpha_Kca31 -
    Ybeta_Kca31*Y_Kca31)*allow_change;
  states[36] = dt*dY_Kca31_dt + Y_Kca31;
  const double I_Kca31 = g_Kca31*(-E_K + v)*Y_Kca31;

  // Expressions for the Cav2.1 component
  const double q10_Cav21 = 3.;
  const double qt_Cav21 = std::pow(q10_Cav21, -2.3 + 0.1*Temp);
  const double vhalfm_Cav21 = -29.4580000000000;
  const double cvm_Cav21 = 8.42900000000000;
  const double vshift_Cav21 = 0.;
  const double minf_Cav21 = 1.0/(1.0 + std::exp((vhalfm_Cav21 + vshift_Cav21 -
    v)/cvm_Cav21));
  const double taum_Cav21 = 1.0*(v >= -40. ? 0.2702 +
    1.1622*std::exp(0.00609050490285645*(26.798 + v)*(-26.798 - v)) :
    0.6923*std::exp(0.00091796007240869*v))/qt_Cav21;
  const double dm_Cav21_dt = (-m_Cav21 + minf_Cav21)*allow_change/taum_Cav21;
  states[37] = dt*dm_Cav21_dt + m_Cav21;
  const double z_Ca = 2.;
  const double zeta_Ca = FARADAY*z_Ca*v/(R*(273.19 + Temp));
  const double ghk = (std::fabs(1. - std::exp(-zeta_Ca)) < 1.0e-6 ?
    1.0e-6*F*z_Ca*(1. + zeta_Ca/2.)*(cai - cao*std::exp(-zeta_Ca)) :
    1.0e-6*F*z_Ca*(cai - cao*std::exp(-zeta_Ca))*zeta_Ca/(1. -
    std::exp(-zeta_Ca)));
  const double I_Cav21 = 1000.0*g_Cav21*(m_Cav21*m_Cav21*m_Cav21)*ghk;

  // Expressions for the Cav3.1 component
  const double q10_Cav31 = 3.;
  const double qt_Cav31 = std::pow(q10_Cav31, -3.7 + 0.1*Temp);
  const double v0_m_inf_Cav31 = -52.0000000000000;
  const double v0_h_inf_Cav31 = -72.0000000000000;
  const double k_m_inf_Cav31 = -5.00000000000000;
  const double k_h_inf_Cav31 = 7.00000000000000;
  const double C_tau_m_Cav31 = 1.00000000000000;
  const double A_tau_m_Cav31 = 1.00000000000000;
  const double v0_tau_m1_Cav31 = -40.0000000000000;
  const double v0_tau_m2_Cav31 = -102.000000000000;
  const double k_tau_m1_Cav31 = 9.00000000000000;
  const double k_tau_m2_Cav31 = -18.0000000000000;
  const double C_tau_h_Cav31 = 15.0000000000000;
  const double A_tau_h_Cav31 = 1.00000000000000;
  const double v0_tau_h1_Cav31 = -32.0000000000000;
  const double k_tau_h1_Cav31 = 7.00000000000000;
  const double minf_Cav31 = 1.0/(1. + std::exp((-v0_m_inf_Cav31 +
    v)/k_m_inf_Cav31));
  const double hinf_Cav31 = 1.0/(1. + std::exp((-v0_h_inf_Cav31 +
    v)/k_h_inf_Cav31));
  const double taum_Cav31 = (v <= -90. ? 1. : (C_tau_m_Cav31 +
    A_tau_m_Cav31/(std::exp((-v0_tau_m1_Cav31 + v)/k_tau_m1_Cav31) +
    std::exp((-v0_tau_m2_Cav31 + v)/k_tau_m2_Cav31)))/qt_Cav31);
  const double tauh_Cav31 = (C_tau_h_Cav31 +
    A_tau_h_Cav31*std::exp(-(-v0_tau_h1_Cav31 + v)/k_tau_h1_Cav31))/qt_Cav31;
  const double dm_Cav31_dt = (-m_Cav31 + minf_Cav31)*allow_change/taum_Cav31;
  states[38] = dt*dm_Cav31_dt + m_Cav31;
  const double dh_Cav31_dt = (-h_Cav31 + hinf_Cav31)*allow_change/tauh_Cav31;
  states[39] = dt*dh_Cav31_dt + h_Cav31;
  const double I_Cav31 = 1000.0*g_Cav31*(m_Cav31*m_Cav31)*ghk*h_Cav31;

  // Expressions for the Cav3.2 component
  const double phi_m_Cav32 = 6.89864830730607;
  const double phi_h_Cav32 = 3.73719281884655;
  const double m_inf_Cav32 = 1.0/(1. +
    0.000607957605998435*std::exp(-0.135135135135135*v));
  const double h_inf_Cav32 = 1.0/(1. +
    148461.062050627*std::exp(0.139275766016713*v));
  const double tau_m_Cav32 = (1.9 +
    1.0/(22.404093719072*std::exp(0.0840336134453781*v) +
    0.00189854653589182*std::exp(-v/21.)))/phi_m_Cav32;
  const double tau_h_Cav32 = 13.7 + (1942.0 +
    55178666.3027343*std::exp(0.108695652173913*v))/(phi_h_Cav32*(1. +
    30321871936.1927*std::exp(0.27027027027027*v)));
  const double dm_Cav32_dt = (-m_Cav32 + m_inf_Cav32)*allow_change/tau_m_Cav32;
  states[40] = dt*dm_Cav32_dt + m_Cav32;
  const double dh_Cav32_dt = (-h_Cav32 + h_inf_Cav32)*allow_change/tau_h_Cav32;
  states[41] = dt*dh_Cav32_dt + h_Cav32;
  const double E_Ca = R*(273.15 + Temp)*std::log(cao/cai)/(2.*FARADAY);
  const double I_Cav32 = g_Cav32*(m_Cav32*m_Cav32)*(-E_Ca + v)*h_Cav32;

  // Expressions for the Cav3.3 component
  const double q10_Cav33 = 2.30000000000000;
  const double qt_Cav33 = std::pow(q10_Cav33, -2.8 + 0.1*Temp);
  const double gCav3_3bar = 1.00000000000000e-5;
  const double vhalfn_Cav33 = -41.5000000000000;
  const double vhalfl_Cav33 = -69.8000000000000;
  const double kn_Cav33 = 6.20000000000000;
  const double kl_Cav33 = -6.10000000000000;
  const double n_inf_Cav33 = 1.0/(1.0 + std::exp((vhalfn_Cav33 - v)/kn_Cav33));
  const double l_inf_Cav33 = 1.0/(1.0 + std::exp((vhalfl_Cav33 - v)/kl_Cav33));
  const double tau_n_Cav33 = (v > -60. ? 7.2 +
    0.02*std::exp(-0.0680272108843537*v) : 79.5 +
    2.0*std::exp(-0.10752688172043*v))/qt_Cav33;
  const double tau_l_Cav33 = (v > -60. ? 0.875*std::exp(120./41. + v/41.) :
    260.0)/qt_Cav33;
  const double w_Cav33 = FARADAY*z_Ca*v/(R*(273.14 + Temp));
  const double ghk_Cav33 = -FARADAY*z_Ca*(cao -
    cai*std::exp(w_Cav33))*w_Cav33/(-1.0 + std::exp(w_Cav33));
  const double dn_Cav33_dt = (-n_Cav33 + n_inf_Cav33)*allow_change/tau_n_Cav33;
  states[42] = dt*dn_Cav33_dt + n_Cav33;
  const double dl_Cav33_dt = (-l_Cav33 + l_inf_Cav33)*allow_change/tau_l_Cav33;
  states[43] = dt*dl_Cav33_dt + l_Cav33;
  const double I_Cav33 =
    gCav3_3bar*g_Cav33*(n_Cav33*n_Cav33*n_Cav33)*ghk_Cav33*l_Cav33;

  // Expressions for the HCN1 component
  const double q10_HCN1 = 3.00000000000000;
  const double qt_HCN1 = std::pow(q10_HCN1, -3.7 + 0.1*Temp);
  const double ratetau_HCN1 = 1.00000000000000;
  const double ljp_HCN1 = 9.30000000000000;
  const double v_inf_half_noljp_HCN1 = -90.3000000000000;
  const double v_inf_k_HCN1 = 9.67000000000000;
  const double v_tau_const_HCN1 = 0.00180000000000000;
  const double v_tau_half1_noljp_HCN1 = -68.0000000000000;
  const double v_tau_half2_noljp_HCN1 = -68.0000000000000;
  const double v_tau_k1_HCN1 = -22.0000000000000;
  const double v_tau_k2_HCN1 = 7.14000000000000;
  const double v_inf_half_HCN1 = v_inf_half_noljp_HCN1 - ljp_HCN1;
  const double v_tau_half1_HCN1 = v_tau_half1_noljp_HCN1 - ljp_HCN1;
  const double v_tau_half2_HCN1 = v_tau_half2_noljp_HCN1 - ljp_HCN1;
  const double hinf_HCN1 = 1.0/(1. + std::exp((-v_inf_half_HCN1 +
    v)/v_inf_k_HCN1));
  const double tauh_HCN1 =
    ratetau_HCN1/(v_tau_const_HCN1*(std::exp((-v_tau_half1_HCN1 +
    v)/v_tau_k1_HCN1) + std::exp((-v_tau_half2_HCN1 +
    v)/v_tau_k2_HCN1))*qt_HCN1);
  const double dh_HCN1_dt = (-h_HCN1 + hinf_HCN1)*allow_change/tauh_HCN1;
  states[44] = dt*dh_HCN1_dt + h_HCN1;
  const double I_HCN1 = g_HCN1*(-E_h + v)*h_HCN1;

  // Expressions for the leak component
  const double I_leak = g_leak*(-E_leak + v);

  // Expressions for the Membrane component
  const double I_tot = I_Cav21 + I_Cav31 + I_Cav32 + I_Cav33 + I_HCN1 +
    I_Kca11 + I_Kca22 + I_Kca31 + I_Kir2x + I_Kv11 + I_Kv15 + I_Kv33 + I_Kv34 +
    I_Kv43 + I_Na + I_leak;
  const double dv_dt = -1.0*I_tot*allow_change/Cm;
  states[45] = dt*dv_dt + v;
}
