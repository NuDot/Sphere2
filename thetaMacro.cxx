{

#include <math.h>
double px[20000], py[20000], pz[20000], track[20000];
  epgTree->SetBranchAddress("epg_pxi",&px);
  epgTree->SetBranchAddress("epg_pyi",&py);
  epgTree->SetBranchAddress("epg_pzi",&pz);
  epgTree->SetBranchAddress("epg_tIDi",&track); 

 
for( int n_entry = 0; n_entry < 100; n_entry++ ) { 
  epgTree->GetEntry(n_entry);
  std::cout<<"---------------------EVT NUM: "<<n_entry<<std::endl;  
  TVector3 vec1, vec2; 
  for( int particle = 0; particle < 20000; particle++ ) {
    std::cout<<"TRACK ID: "<<track[particle]<<std::endl; 
    if( track[particle] == 1 ) { 
      vec1 = TVector3(px[particle], py[particle], pz[particle]); 
  //    std::cout<<"FOUND TRACK ID 1"<<std::endl; 
    }
    else if( track[particle] == 2 ) {
      vec2 = TVector3(px[particle], py[particle], pz[particle]);
//      std::cout<<"FOUND TRACK ID 2"<<std::endl; 
    }
  }
}

}

