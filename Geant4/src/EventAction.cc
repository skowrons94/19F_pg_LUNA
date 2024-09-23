/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "Analysis.hh"
#include "GeometryManager.hh"

#include "geometry/BGODetector.hh"

#include "geometry/GeometryObject.hh"

#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4EventManager.hh"
#include "G4PrimaryParticle.hh"

#include "G4ThreeVector.hh"

#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Let GeometryManager loop over (sensitive) geometry objects and let them
    // fill their outputs.
    auto *gm = GeometryManager::GetInstance();
    gm->FillOutputs(event);

    TrackAnnihilation(event);
    
}


void EventAction::TrackAnnihilation(const G4Event* event)
{

    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

    for(G4int i=0; i<n_trajectories; i++) 
    {
      G4Trajectory* trj = (G4Trajectory*)((*(event->GetTrajectoryContainer()))[i]);
      if( trj->GetParticleName() == "gamma" )
	{
	  if( trj->GetInitialKineticEnergy() > 0.51 && trj->GetInitialKineticEnergy() < 0.52 )
	    {
	      G4ThreeVector vec = trj->GetPoint(0)->GetPosition();
	      
	      auto *gm = GeometryManager::GetInstance();
	      GeometryObject *bgo = gm->GetGeometry("BGODetector");
	      bgo->FillPosition(vec);
	      break;
	    }
	}
    }
  
}
