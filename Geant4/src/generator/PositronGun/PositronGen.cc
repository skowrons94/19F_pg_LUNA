/// \file IsotropicGunGen.cc
/// \brief Implementation of the IsotropicGunGen class

#include "generator/PositronGun/PositronGunGen.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"

#include "TGraph.h"
#include "TSpline.h"
#include "TRandom.h"

#include <memory>
using std::make_shared;

#include <stdexcept>
using std::runtime_error;

PositronGunGen::PositronGunGen()
    : G4VUserPrimaryGeneratorAction(),
      G4UImessenger(),
      fParticleGun(nullptr),
      m_position(0, 0, 0)
{

    m_setPositionCmd = make_shared<G4UIcmdWith3VectorAndUnit>("/PrimaryGenerator/PositronGun/position", this);
    m_setPositionCmd->SetGuidance("Set position of primary vertex.");
    m_setPositionCmd->SetParameterName("x", "y", "z", false);
    m_setPositionCmd->SetUnitCategory("Length");

    m_setInputFileNameCmd = make_shared<G4UIcmdWithAString>("/PrimaryGenerator/PositronGun/file", this);
    m_setInputFileNameCmd->SetGuidance("Set the name of the input file with the decay spectrum.");
    m_setInputFileNameCmd->SetParameterName("file name", false);

    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);

    rand = new TRandom( );
    
    fParticleGun->SetParticleDefinition(G4Positron::Positron());

}

PositronGunGen::~PositronGunGen()
{
    delete fParticleGun;
}

void PositronGunGen::GeneratePrimaries(G4Event* anEvent)
{
    // This function is called at the begining of event

    fParticleGun->SetParticlePosition( m_position );

    do
      {
	  energy = rand->Rndm( )*m_energy;
	  prob = rand->Rndm(  );
      } while( prob > m_pdf->Eval( energy ) );
    
    fParticleGun->SetParticleEnergy(energy*CLHEP::MeV);
    // Set gun direction randomly
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

    // Generate the vertex (position + particle) for the event
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PositronGunGen::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == m_setPositionCmd.get())
    {
        m_position = m_setPositionCmd->GetNew3VectorValue(newValue);
    }
    else if(command == m_setInputFileNameCmd.get())
    {
        m_graph = new TGraph(newValue.c_str());
        m_pdf = new TSpline3( "pdf", m_graph );
        m_energy = m_graph->GetPointX( m_graph->GetMaxSize( ) - 1 );
    }
    else
    {
        throw runtime_error("Unknown command in GammaDecaySchemeGen::SetNewValue()");
    }
}