//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef NuclideGunGen_h
#define NuclideGunGen_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "TRandom.h"

#include "G4UImessenger.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using CLHEP::mm;

using std::shared_ptr;

class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADouble;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NuclideGunGen : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
  public:
    NuclideGunGen();    
   ~NuclideGunGen();

  public:
    virtual void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() { return fParticleGun;} ;

    void SetNewValue(G4UIcommand* command, G4String newValue);
            
  private:
    G4ParticleGun*  fParticleGun;
    G4ParticleGun*  fGammaGun;

    G4ThreeVector m_position;
    G4ThreeVector m_position_rand; // position of primary vertex

    G4double beam_sigma = 5 * CLHEP::mm;
    
    TRandom* rand;
    G4double theta;
    G4double offset;
    G4double offsetX;
    G4double offsetY;

    G4double m_mass = 137;
    G4double m_charge = 55;

    G4double m_gamma   = 0;
    G4double m_gammaEnergy = 0;

    shared_ptr<G4UIcmdWithADouble> m_selectMassCmd;
    shared_ptr<G4UIcmdWithADouble> m_selectChargeCmd;
    shared_ptr<G4UIcmdWith3VectorAndUnit> m_setPositionCmd;
    shared_ptr<G4UIcmdWithADouble> m_setGammaCmd;
    shared_ptr<G4UIcmdWithADouble> m_setGammaEnergyCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
