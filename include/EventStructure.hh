#ifndef EventStructure_H
#define EventStructure_H 1

const int MAX_epg=10000000;//1000;
const int MAX_phot=10000000;

struct event
{
  int evt_num;

  float C10dT;
  float B10dT;
  double t0;

  int isB10_1740;
  int isB10_718;
  double C10_0_ti;
  double nue_ti;
  double B10_0_ti;
  double B10_718_ti;
  double B10_1740_ti;

  int N_epg_i;
  int epg_pIDi[MAX_epg];
  int epg_tIDi[MAX_epg];
  float epg_qi[MAX_epg];
  double epg_ti[MAX_epg];
  float epg_xi[MAX_epg];
  float epg_yi[MAX_epg];
  float epg_zi[MAX_epg];
  float epg_pxi[MAX_epg];
  float epg_pyi[MAX_epg];
  float epg_pzi[MAX_epg];
  float epg_ei[MAX_epg];
  int epg_isOPi[MAX_epg];
  int epg_isChei[MAX_epg];

  int N_epg_f;
  int epg_pIDf[MAX_epg];
  int epg_tIDf[MAX_epg];
  float epg_qf[MAX_epg];
  double epg_tf[MAX_epg]; //keep time as double for better precision
  double epg_ltf[MAX_epg];
  float epg_xf[MAX_epg];
  float epg_yf[MAX_epg];
  float epg_zf[MAX_epg];
  float epg_pxf[MAX_epg];
  float epg_pyf[MAX_epg];
  float epg_pzf[MAX_epg];
  float epg_ef[MAX_epg];
  float epg_lf[MAX_epg];
  int epg_isOPf[MAX_epg];
  int epg_isChef[MAX_epg];

  float edep;
  float edep_cor;

  int N_phot;
  float t_start[MAX_phot];
  float x_start[MAX_phot];
  float y_start[MAX_phot];
  float z_start[MAX_phot];
  float x_hit[MAX_phot];
  float y_hit[MAX_phot];
  float z_hit[MAX_phot];
  float x_vtx[MAX_phot];
  float y_vtx[MAX_phot];
  float z_vtx[MAX_phot];
  float cos_theta[MAX_phot];
  float photon_wavelength[MAX_phot];
  float true_time[MAX_phot];
  int PE_creation[MAX_phot];
  float PE_time[MAX_phot];
  int detector_coverage_included[MAX_phot];
  float true_time_corrected[MAX_phot];
  float PE_time_corrected[MAX_phot];
  float cos_theta_reco[MAX_phot];
  float theta_reco[MAX_phot];
  float phi_reco[MAX_phot];
  int process[MAX_phot];  
 
  double trueVtxX;
  double trueVtxY;
  double trueVtxZ; 
  double trueDirX;
  double trueDirY;
  double trueDirZ;

};

#endif
