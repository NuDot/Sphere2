{

#include <math.h>
float px[500000], py[500000], pz[500000];
int track[500000];
epgTree->SetBranchAddress("epg_pxi",&px);
epgTree->SetBranchAddress("epg_pyi",&py);
epgTree->SetBranchAddress("epg_pzi",&pz);
epgTree->SetBranchAddress("epg_tIDi",&track); 

TH1* angles = new TH1F("angles", "angles", 100, 0.00, 3.25);
 
for( int n_entry = 0; n_entry < 1000; n_entry++ ) { 
  epgTree->GetEntry(n_entry);
  std::cout<<"EVT NUM: "<<n_entry<<std::endl;  
  TVector3 vec1, vec2; 
  for( int particle = 0; particle < 500000; particle++ ) {
    if( track[particle] == 1 ) { 
      std::cout<<"Particle 1 Momentum: "<<px[particle]<<"   "<<py[particle]<<"   "<<pz[particle]<<std::endl; 
      vec1 = TVector3(px[particle], py[particle], pz[particle]); 
//      std::cout<<"TRACKING NUMBER 1 FOUND"<<std::endl; 
    }
    else if( track[particle] == 2 ) {
      std::cout<<"Particle 2 Momentum: "<<px[particle]<<"   "<<py[particle]<<"   "<<pz[particle]<<std::endl;
      vec2 = TVector3(px[particle], py[particle], pz[particle]);
//      std::cout<<"TRACKING NUMBER 2 FOUND"<<std::endl;
    }
  }
  std::cout<<"Angle Between Particles 1 and 2: "<<vec1.Angle(vec2)<<std::endl;
  angles->Fill(vec1.Angle(vec2)); 
}

angles->Draw(); 

}

