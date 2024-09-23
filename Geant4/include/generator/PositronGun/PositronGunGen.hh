/// \file PositronGen.hh
/// \brief Definition of the PositronGen class

#ifndef PositronGunGen_h
#define PositronGunGen_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"
#include "globals.hh"

#include "TGraph.h"
#include "TRandom.h"
#include "TSpline.h"

class G4ParticleGun;
class G4Event;

class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;

using std::shared_ptr;

/// The primary generator action class with particle gun.
///
/// At the moment it only instantiates the particle gun and uses it to create a
/// single particle per event.
/// This means that the /gun/ commands can be used in the macro to set the
/// energy and position of the positron.
/// The direction of the particle is sampled randomly for every event, thus
/// setting the direction in the macro using /gun/direction will be disregarded.


class PositronGunGen : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
public:
  PositronGunGen();
  virtual ~PositronGunGen();

  virtual void GeneratePrimaries(G4Event* event);

  void SetNewValue(G4UIcommand* command, G4String newValue);
  
private:
  G4ParticleGun*  fParticleGun;

  G4ThreeVector m_position;
  G4double m_energy;
  G4String m_file;
  TGraph* m_graph;
  TSpline3* m_pdf;
  
  TRandom* rand;

  G4double energy;
  G4double prob;

  shared_ptr<G4UIcmdWith3VectorAndUnit> m_setPositionCmd;
  shared_ptr<G4UIcmdWithAString>        m_setInputFileNameCmd;

};

#endif
