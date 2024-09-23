/// \file IsotropicGunGen.cc
/// \brief Implementation of the IsotropicGunGen class

#include "generator/IsotropicGun/IsotropicGunGen.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <math.h>

using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;

#include <memory>
using std::make_shared;

#include <stdexcept>
using std::runtime_error;

IsotropicGunGen::IsotropicGunGen()
    : G4VUserPrimaryGeneratorAction(),
      G4UImessenger(), 
      fParticleGun(nullptr), 
      m_position(0, 0, 0)
{
    m_setPositionCmd = make_shared<G4UIcmdWith3VectorAndUnit>("/PrimaryGenerator/IsotropicGun/position", this);
    m_setPositionCmd->SetGuidance("Set position of primary vertex.");
    m_setPositionCmd->SetParameterName("x", "y", "z", false);
    m_setPositionCmd->SetUnitCategory("Length");

    m_selectEnergyCmd = make_shared<G4UIcmdWithADoubleAndUnit>("/PrimaryGenerator/IsotropicGun/energy", this);
    m_selectEnergyCmd->SetGuidance("Select the energy of the gamma.");
    m_selectEnergyCmd->SetParameterName("energy", false);
    m_selectEnergyCmd->SetUnitCategory("Energy");

    m_selectNParticlesCmd = make_shared<G4UIcmdWithADouble>("/PrimaryGenerator/IsotropicGun/number", this);
    m_selectNParticlesCmd->SetGuidance("Select the number of the gammas.");
    m_selectNParticlesCmd->SetParameterName("N", false);
    
}

IsotropicGunGen::~IsotropicGunGen()
{
    delete fParticleGun;
}

void IsotropicGunGen::GeneratePrimaries(G4Event* anEvent)
{
    // This function is called at the begining of event

    fParticleGun->SetParticleDefinition(G4Gamma::Gamma());

    fParticleGun->SetParticleEnergy(m_energy*CLHEP::MeV);
    fParticleGun->SetParticlePosition(m_position);

    // Set gun direction randomly
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

    // Generate the vertex (position + particle) for the event
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

void IsotropicGunGen::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == m_setPositionCmd.get())
    {
        m_position = m_setPositionCmd->GetNew3VectorValue(newValue);
    }
    else if(command == m_selectEnergyCmd.get())
    {
        m_energy = m_selectEnergyCmd->GetNewDoubleValue(newValue);
    }
    else if(command == m_selectNParticlesCmd.get())
    {
        m_number = m_selectNParticlesCmd->GetNewDoubleValue(newValue);
	fParticleGun = new G4ParticleGun(m_number);
    }
    else
    {
        throw runtime_error("Unknown command in GammaDecaySchemeGen::SetNewValue()");
    }
}
