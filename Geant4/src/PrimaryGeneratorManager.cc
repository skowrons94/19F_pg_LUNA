/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorManager.hh"

#include "generator/IsotropicGun/IsotropicGunGen.hh"
#include "generator/GammaDecayScheme/GammaDecaySchemeGen.hh"
#include "generator/PositronGun/PositronGunGen.hh"
#include "generator/NuclideGun/NuclideGunGen.hh"

#include "G4Event.hh"
#include "G4UIcmdWithAString.hh"

#include <stdexcept>
using std::runtime_error;

PrimaryGeneratorManager::PrimaryGeneratorManager()
    : G4VUserPrimaryGeneratorAction(),
      G4UImessenger()
{
    /// Define macro command
    m_selectPGcmd = make_shared<G4UIcmdWithAString>("/PrimaryGenerator/select", this);
    m_selectPGcmd->SetGuidance("Choose primary generator.");
    m_selectPGcmd->SetParameterName("Primary generator name.", false);
    m_selectPGcmd->SetCandidates("IsotropicGun GammaDecayScheme PositronGun NuclideGun");

    /// Initialize primary generators
    m_pgIsotropicGun     = make_shared<IsotropicGunGen>();
    m_pgGammaDecayScheme = make_shared<GammaDecaySchemeGen>();
    m_pgPositronGun = make_shared<PositronGunGen>();
    m_pgNuclideGun = make_shared<NuclideGunGen>();
}


void PrimaryGeneratorManager::GeneratePrimaries(G4Event* anEvent)
{
    switch (m_selectedPG)
    {
        case pgIsotropicGun:
            m_pgIsotropicGun->GeneratePrimaries(anEvent);
            break;

        case pgGammaDecayScheme:
            m_pgGammaDecayScheme->GeneratePrimaries(anEvent);
            break;

        case pgPositronGun:
            m_pgPositronGun->GeneratePrimaries(anEvent);
            break;

	case pgNuclideGun:
            m_pgNuclideGun->GeneratePrimaries(anEvent);
            break;

        case pgUndefined:
            throw runtime_error("No primary generator selected!");
            break;

        default:
            throw runtime_error("Unhandled case in PrimaryGeneratorManager::GeneratePrimaries().");
            break;
    }
}

void PrimaryGeneratorManager::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == m_selectPGcmd.get())
    {
        if (newValue == "IsotropicGun")
        {
            m_selectedPG = pgIsotropicGun;
        }
        else if (newValue == "GammaDecayScheme")
        {
            m_selectedPG = pgGammaDecayScheme;
        }
        else if (newValue == "PositronGun")
        {
            m_selectedPG = pgPositronGun;
        }
	else if (newValue == "NuclideGun")
        {
            m_selectedPG = pgNuclideGun;
        }
        else
        {
            G4cerr << "Unknown primary generator to be selected: " << newValue << G4endl;
            throw runtime_error("Unknown primary generator in PrimaryGeneratorManager::SetNewValue().");
        }
    }
    else
    {
        throw runtime_error("Unhandled command in PrimaryGeneratorManager::SetNewValue().");
    }
}

