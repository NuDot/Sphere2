#include "Sphere1TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"

Sphere1TrackingAction::Sphere1TrackingAction(event* fEv)
{
//  tTree = fTree1;
//  tTree2 = fTree2;
   tEv = fEv;
}
 
void Sphere1TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="opticalphoton"&&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="gamma"&&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="e-"&&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="e+" &&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="C10[0.0]" &&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="C11[0.0]" &&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="Tl208[0.0]" &&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="U238[0.0]" &&
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())!="Th232[0.0]")
  {
    G4cout<<"NAME = "<<aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName()<<G4endl;
    G4cout<<"Proc_name = "<<aTrack->GetCreatorProcess()->GetProcessName()<<G4endl;
    G4cout<<"vtx_pos: r = "<<sqrt(aTrack->GetVertexPosition().getX()*aTrack->GetVertexPosition().getX()+
                                  aTrack->GetVertexPosition().getY()*aTrack->GetVertexPosition().getY()+
				  aTrack->GetVertexPosition().getZ()*aTrack->GetVertexPosition().getZ() )<<G4endl;
  }

  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="C10[0.0]")
  {
    tEv->C10_0_ti=aTrack->GetGlobalTime();
    int prec=G4cout.precision();
    G4cout.precision(15);
    G4cout<<"time_C10_0 = "<<tEv->C10_0_ti<<G4endl;
    G4cout.precision(prec);
  }

  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="nu_e")
  {
    tEv->nue_ti=aTrack->GetGlobalTime();
//    if(aTrack->GetDynamicParticle()->GetKineticEnergy()<0.6) 
    tEv->t0=tEv->nue_ti;
    int prec=G4cout.precision();
    G4cout.precision(15);
    G4cout<<"time_nue = "<<tEv->nue_ti<<G4endl;
    G4cout.precision(prec);
    tEv->C10dT=tEv->nue_ti-tEv->C10_0_ti;
    G4Track* bTrack;
    bTrack = (G4Track*)aTrack;
    bTrack->SetTrackStatus(fStopAndKill);
  }

  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="anti_nu_e")
  {
    tEv->t0=aTrack->GetGlobalTime();
    G4Track* bTrack;
    bTrack = (G4Track*)aTrack;
    bTrack->SetTrackStatus(fStopAndKill);
  }


  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="B10[718.3]" || 
     (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="B10[718.5]" )
  {
    tEv->isB10_718=1;
    tEv->B10_718_ti=aTrack->GetGlobalTime();
    int prec=G4cout.precision();
    G4cout.precision(15);
    G4cout<<"time_B10_718 = "<<tEv->B10_718_ti<<G4endl;
    G4cout.precision(prec);
  }
  
  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="B10[1740.2]")
  {
    tEv->isB10_1740=1;
    tEv->B10_1740_ti=aTrack->GetGlobalTime();
    int prec=G4cout.precision();
    G4cout.precision(15);
    G4cout<<"timee_B10_1740 = "<<tEv->B10_1740_ti<<G4endl;
    G4cout.precision(prec);
  }

  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="B10[0.0]")
  {
    tEv->B10_0_ti=aTrack->GetGlobalTime();
    int prec=G4cout.precision();
    G4cout.precision(15);
    G4cout<<"time_B10_0 = "<<tEv->B10_0_ti<<G4endl;
    G4cout<<"B10dT = "<<-(tEv->nue_ti-tEv->B10_0_ti)<<G4endl;
    G4cout.precision(prec);
    tEv->B10dT=-(tEv->nue_ti-tEv->B10_0_ti);
  }
/*
  if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="opticalphoton")
  {

    G4Track* bTrack;
    bTrack = (G4Track*)aTrack;
    bTrack->SetTrackStatus(fStopAndKill);
  }
*/  
/*	if((aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName())=="opticalphoton")
	{
	//	    tEv->NPhotons++;
    tEv->Energy[tEv->NPhotons]=aTrack->GetDynamicParticle()->GetKineticEnergy();
    tEv->X[tEv->NPhotons]=aTrack->GetPosition().x();
    tEv->Y[tEv->NPhotons]=aTrack->GetPosition().y();
    tEv->Z[tEv->NPhotons]=aTrack->GetPosition().z();
    tEv->Px[tEv->NPhotons]=aTrack->GetMomentumDirection().x();
    tEv->Py[tEv->NPhotons]=aTrack->GetMomentumDirection().y();
    tEv->Pz[tEv->NPhotons]=aTrack->GetMomentumDirection().z();
    tEv->Tg[tEv->NPhotons]=aTrack->GetGlobalTime();
    tEv->Tl[tEv->NPhotons]=aTrack->GetLocalTime();
    tEv->ParentID[tEv->NPhotons] = aTrack->GetParentID();
    G4cout<<"ParentID = "<<aTrack->GetParentID()<<G4endl;
		tEv->NPhotons++;
//		G4cout << aTrack->GetPosition().x() << G4endl;
//		G4Track* bTrack;
//		bTrack = (G4Track*)aTrack;
//		bTrack->SetTrackStatus(fStopAndKill);
	}
*/
  G4String pName = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if(pName=="e+"||pName=="e-"||pName=="gamma"||pName=="opticalphoton")
  {
    int id=aTrack->GetTrackID();
    bool new_track=true;
    for(int i=0;i!=tEv->N_epg_i;++i)
    {
      if(id==tEv->epg_tIDi[i]) new_track=false;
    }
    if(new_track)
    {
      tEv->epg_tIDi[tEv->N_epg_i]=aTrack->GetTrackID();
      tEv->epg_pIDi[tEv->N_epg_i]=aTrack->GetParentID();
      tEv->epg_qi[tEv->N_epg_i]=aTrack->GetDynamicParticle()->GetCharge();
      tEv->epg_ti[tEv->N_epg_i]=aTrack->GetGlobalTime();
      tEv->epg_xi[tEv->N_epg_i]=aTrack->GetPosition().x();
      tEv->epg_yi[tEv->N_epg_i]=aTrack->GetPosition().y();
      tEv->epg_zi[tEv->N_epg_i]=aTrack->GetPosition().z();
      tEv->epg_pxi[tEv->N_epg_i]=aTrack->GetMomentumDirection().x();
      tEv->epg_pyi[tEv->N_epg_i]=aTrack->GetMomentumDirection().y();
      tEv->epg_pzi[tEv->N_epg_i]=aTrack->GetMomentumDirection().z();
      tEv->epg_ei[tEv->N_epg_i]=aTrack->GetDynamicParticle()->GetKineticEnergy();
      tEv->epg_isOPi[tEv->N_epg_i]=0;
      tEv->epg_isChei[tEv->N_epg_i]=0;

      if(pName=="opticalphoton") 
      {
        tEv->epg_isOPi[tEv->N_epg_i]=1;
        if(aTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov")
  	  tEv->epg_isChei[tEv->N_epg_i]=1;

//	G4Track* bTrack;
//	bTrack = (G4Track*)aTrack;
//	bTrack->SetTrackStatus(fStopAndKill);
      }

      tEv->N_epg_i++;    
    }
  }
/* AE (09/10/2014) I don't remember why I had to separate photons this way...
  if(pName=="opticalphoton")
  {
      tEv->epg_isOP[tEv->N_epg_i]=1;
      if(aTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov")
      {
	tEv->epg_isChe[tEv->N_epg_i]=1;
      } else {
	tEv->epg_isChe[tEv->N_epg_i]=0;
      }

      tEv->epg_tIDi[tEv->N_epg_i]=aTrack->GetTrackID();
      tEv->epg_pIDi[tEv->N_epg_i]=aTrack->GetParentID();
      tEv->epg_qi[tEv->N_epg_i]=aTrack->GetDynamicParticle()->GetCharge();
      tEv->epg_ti[tEv->N_epg_i]=aTrack->GetGlobalTime();
      tEv->epg_xi[tEv->N_epg_i]=aTrack->GetPosition().x();
      tEv->epg_yi[tEv->N_epg_i]=aTrack->GetPosition().y();
      tEv->epg_zi[tEv->N_epg_i]=aTrack->GetPosition().z();
      tEv->epg_pxi[tEv->N_epg_i]=aTrack->GetMomentumDirection().x();
      tEv->epg_pyi[tEv->N_epg_i]=aTrack->GetMomentumDirection().y();
      tEv->epg_pzi[tEv->N_epg_i]=aTrack->GetMomentumDirection().z();
      tEv->epg_ei[tEv->N_epg_i]=aTrack->GetDynamicParticle()->GetKineticEnergy();

    tEv->N_epg_i++;
//    G4Track* bTrack;
//    bTrack = (G4Track*)aTrack;
//    bTrack->SetTrackStatus(fStopAndKill);
  }
*/

  return;
}  

void Sphere1TrackingAction::PostUserTrackingAction(const G4Track* aTrack) 
{
  G4String pName = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if(pName=="e+"||pName=="e-"||pName=="gamma"||pName=="opticalphoton")
  {
    int id=aTrack->GetTrackID();
    bool new_track=true;
    int n=0;
    for(int i=0;i!=tEv->N_epg_f;++i)
    {
      if(id==tEv->epg_tIDf[i]) 
      { 
	new_track=false;
	n=i;
      }
    }
    if(new_track)
    {
      tEv->epg_tIDf[tEv->N_epg_f]=aTrack->GetTrackID();
      tEv->epg_pIDf[tEv->N_epg_f]=aTrack->GetParentID();
      tEv->epg_qf[tEv->N_epg_f]=aTrack->GetDynamicParticle()->GetCharge();
      tEv->epg_tf[tEv->N_epg_f]=aTrack->GetGlobalTime();
      tEv->epg_ltf[tEv->N_epg_f]=aTrack->GetLocalTime();
      tEv->epg_xf[tEv->N_epg_f]=aTrack->GetPosition().x();
      tEv->epg_yf[tEv->N_epg_f]=aTrack->GetPosition().y();
      tEv->epg_zf[tEv->N_epg_f]=aTrack->GetPosition().z();
      tEv->epg_pxf[tEv->N_epg_f]=aTrack->GetMomentumDirection().x();
      tEv->epg_pyf[tEv->N_epg_f]=aTrack->GetMomentumDirection().y();
      tEv->epg_pzf[tEv->N_epg_f]=aTrack->GetMomentumDirection().z();
      tEv->epg_ef[tEv->N_epg_f]=aTrack->GetDynamicParticle()->GetKineticEnergy();
      tEv->epg_lf[tEv->N_epg_f]=aTrack->GetTrackLength();

      tEv->epg_isOPf[tEv->N_epg_f]=0;
      tEv->epg_isChef[tEv->N_epg_f]=0;

      if(pName=="opticalphoton")
      {
        tEv->epg_isOPf[tEv->N_epg_f]=1;
        if(aTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov")
          tEv->epg_isChef[tEv->N_epg_f]=1;
      }

      tEv->N_epg_f++;

    } else
    {
      if(tEv->epg_tIDf[n]!=aTrack->GetTrackID() ||
         tEv->epg_pIDf[n]!=aTrack->GetParentID() ||
         tEv->epg_qf[n]!=aTrack->GetDynamicParticle()->GetCharge())
        G4cout<<"Particle Mismatch Check TrackingAction!!!"<<G4endl;
      tEv->epg_xf[n]=aTrack->GetPosition().x();
      tEv->epg_yf[n]=aTrack->GetPosition().y();
      tEv->epg_zf[n]=aTrack->GetPosition().z();
      tEv->epg_pxf[n]=aTrack->GetMomentumDirection().x();
      tEv->epg_pyf[n]=aTrack->GetMomentumDirection().y();
      tEv->epg_pzf[n]=aTrack->GetMomentumDirection().z();
      tEv->epg_ef[n]=aTrack->GetDynamicParticle()->GetKineticEnergy();
    }
  }
}
