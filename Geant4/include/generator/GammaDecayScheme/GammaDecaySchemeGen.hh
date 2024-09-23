// GammaDecaySchemeGen.hh

#ifndef GammaDecaySchemeGen_hh
#define GammaDecaySchemeGen_hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"

#include "TRandom.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include <memory>
using std::shared_ptr;

class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class LevelScheme;

using CLHEP::mm;

/*
 * *  Function to read decay information from text file
 *  Format is as follows:
 *  State ENERGY
 *  DAUGHTER    DAUGHTER    DAUGHTER    ...
 *  PROB        PROB        PROB        ...
 *
 *  State NEXT_ENERGY
 *  DAUGHTER    ...
 *  PROB        ...
 *
 *  Lines starting with "//" are skipped
 *  Also everything after "State energy" is ignored
 *  Energies are in keV
 *  Probabilities can be relative or absolute, they will be normalized later
 */


class GammaDecaySchemeGen : public G4VUserPrimaryGeneratorAction, public G4UImessenger
{
public:
    GammaDecaySchemeGen();
    virtual ~GammaDecaySchemeGen() {};

    void GeneratePrimaries(G4Event* anEvent);

    void SetNewValue(G4UIcommand* command, G4String newValue);

private:
    shared_ptr<LevelScheme> m_levels; // Level scheme with selected initial level, used to generate decay gammas

    G4ThreeVector m_position; // position of primary vertex
    G4ThreeVector m_position_rand; // position of primary vertex

    G4double beam_sigma = 6 * mm;

    TRandom* rand;
    G4double theta;
    G4double offset;
    G4double offsetX;
    G4double offsetY;

    shared_ptr<G4UIcmdWith3VectorAndUnit> m_setPositionCmd;
    shared_ptr<G4UIcmdWithAString>        m_setInputFileNameCmd;
    shared_ptr<G4UIcmdWithADoubleAndUnit> m_selectExcitedStateCmd;
};

#endif
