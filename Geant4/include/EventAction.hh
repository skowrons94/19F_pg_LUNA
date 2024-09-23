/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "G4THitsMap.hh"
#include "globals.hh"

/// Event action class

class EventAction : public G4UserEventAction
{
public:
    EventAction() : G4UserEventAction() {}
    virtual ~EventAction() {}

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);

    void TrackAnnihilation(const G4Event* event );
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


