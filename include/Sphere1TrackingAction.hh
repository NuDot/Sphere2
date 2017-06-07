#ifndef Sphere1TrackingAction_H
#define Sphere1TrackingAction_H 1

#include "globals.hh" 
#include "G4UserTrackingAction.hh"
//#include "TTree.h"
#include "EventStructure.hh"

#include <string>

class Sphere1TrackingActionMessenger; 

class Sphere1TrackingAction : public G4UserTrackingAction
{
public:
	Sphere1TrackingAction(event*);
        ~Sphere1TrackingAction();

	void PreUserTrackingAction(const G4Track* aTrack);
	void PostUserTrackingAction(const G4Track* aTrack);

//	double time_C10_0;
//	double time_B10_718;
//	double time_B10_0;
//	double time_nue;
public: 
        //for the Messenger
        void SetKillIon(G4String); 
private:
//	TTree* tTree;
	
        event* tEv;
        
        Sphere1TrackingActionMessenger* tMessenger; 
};
	       
#endif
