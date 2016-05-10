#ifndef Sphere1TrackingAction_H
#define Sphere1TrackingAction_H 1
 
#include "G4UserTrackingAction.hh"
//#include "TTree.h"
#include "EventStructure.hh"

class Sphere1TrackingAction : public G4UserTrackingAction
{
public:
	Sphere1TrackingAction(event*);
        ~Sphere1TrackingAction() {};

	void PreUserTrackingAction(const G4Track* aTrack);
	void PostUserTrackingAction(const G4Track* aTrack);

//	double time_C10_0;
//	double time_B10_718;
//	double time_B10_0;
//	double time_nue;
private:
//	TTree* tTree;
	event* tEv;
};
	       
#endif
