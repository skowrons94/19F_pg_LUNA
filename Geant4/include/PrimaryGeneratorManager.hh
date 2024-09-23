/// \file PrimaryGeneratorManager.hh
/// \brief Definition of the PrimaryGeneratorManager class

#ifndef PrimaryGeneratorManager_h
#define PrimaryGeneratorManager_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>
using std::shared_ptr;
using std::make_shared;

class G4Event;
class G4UIcmdWithAString;

class IsotropicGunGen;
class GammaDecaySchemeGen;
class PositronGunGen;
class NuclideGunGen;

/// The primary generator action manager.
///
/// Allows to select among the primary generator actions to be used.
///


class PrimaryGeneratorManager : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
public:
    PrimaryGeneratorManager();

    virtual ~PrimaryGeneratorManager() {}

    virtual void GeneratePrimaries(G4Event* event);

    void SetNewValue(G4UIcommand* command, G4String newValue);

private:
    shared_ptr<G4UIcmdWithAString> m_selectPGcmd;

  enum {pgIsotropicGun, pgGammaDecayScheme, pgPositronGun, pgNuclideGun, pgUndefined} m_selectedPG = pgUndefined;

    shared_ptr<IsotropicGunGen>     m_pgIsotropicGun;
    shared_ptr<GammaDecaySchemeGen> m_pgGammaDecayScheme;
    shared_ptr<PositronGunGen> m_pgPositronGun;
    shared_ptr<NuclideGunGen> m_pgNuclideGun;

};

#endif
