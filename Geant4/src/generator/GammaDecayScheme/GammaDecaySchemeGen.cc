// GammaDecaySchemeGen.cc

#include "generator/GammaDecayScheme/GammaDecaySchemeGen.hh"

#include "G4Gamma.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RandomDirection.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4Event.hh"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <math.h>

#include "CLHEP/Units/SystemOfUnits.h"

using std::vector;
using std::string;
using std::ifstream;
using std::stringstream;

using CLHEP::keV;

#include "generator/GammaDecayScheme/LevelScheme.hh"

#include <memory>
using std::make_shared;

#include <stdexcept>
using std::runtime_error;


GammaDecaySchemeGen::GammaDecaySchemeGen() : G4VUserPrimaryGeneratorAction(), G4UImessenger(), m_position(0, 0, 0)
{
    m_setPositionCmd = make_shared<G4UIcmdWith3VectorAndUnit>("/PrimaryGenerator/GammaDecayScheme/position", this);
    m_setPositionCmd->SetGuidance("Set position of primary vertex.");
    m_setPositionCmd->SetParameterName("x", "y", "z", false);
    m_setPositionCmd->SetUnitCategory("Length");

    m_setInputFileNameCmd = make_shared<G4UIcmdWithAString>("/PrimaryGenerator/GammaDecayScheme/levelFile", this);
    m_setInputFileNameCmd->SetGuidance("Set the name of the input file with the decay levels.");
    m_setInputFileNameCmd->SetParameterName("file name", false);

    m_selectExcitedStateCmd = make_shared<G4UIcmdWithADoubleAndUnit>("/PrimaryGenerator/GammaDecayScheme/excitedState", this);
    m_selectExcitedStateCmd->SetGuidance("Select the entry level of the scheme.");
    m_selectExcitedStateCmd->SetParameterName("energy", false);
    m_selectExcitedStateCmd->SetUnitCategory("Energy");

    rand = new TRandom( );
}


void GammaDecaySchemeGen::GeneratePrimaries(G4Event* anEvent)
{

    // Sampling the beamspot
    theta = rand->Rndm( )*M_PI;
    offset = rand->Rndm( )*beam_sigma - beam_sigma/2;
    offsetX = offset*std::cos( theta ) * mm;
    offsetY = offset*std::sin( theta ) * mm;

    m_position_rand.setX( m_position.getX( ) + offsetX );
    m_position_rand.setY( m_position.getX( ) + offsetY );
    m_position_rand.setZ( m_position.getZ( ) );
    
    auto *primaryVertex = new G4PrimaryVertex( m_position_rand, 0 ); // t = 0.0

    m_levels->Start();
    while (!m_levels->IsAtEndState())
    {
        auto *primaryParticle = new G4PrimaryParticle(G4Gamma::GammaDefinition());

        primaryParticle->SetMomentumDirection(G4RandomDirection());
        primaryParticle->SetKineticEnergy(m_levels->Decay());

        primaryVertex->SetPrimary(primaryParticle);
    }

    anEvent->AddPrimaryVertex(primaryVertex);
}

void GammaDecaySchemeGen::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == m_setPositionCmd.get())
    {
        m_position = m_setPositionCmd->GetNew3VectorValue(newValue);
    }
    else if(command == m_setInputFileNameCmd.get())
    {
        m_levels = make_shared<LevelScheme>(newValue);
    }
    else if (command == m_selectExcitedStateCmd.get())
    {
        if (m_levels == nullptr)
        {
            throw runtime_error("GammaDecaySchemeGen::SetNewValue(): Tried to select excited state before loading decay scheme.");
        }
        m_levels->SetStartLevel(m_selectExcitedStateCmd->GetNewDoubleValue(newValue));
    }
    else
    {
        throw runtime_error("Unknown command in GammaDecaySchemeGen::SetNewValue()");
    }
}


