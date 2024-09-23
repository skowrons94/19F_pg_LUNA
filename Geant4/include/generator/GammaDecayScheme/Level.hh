// Level.hh

#ifndef Level_hh
#define Level_hh

#include "globals.hh"
#include "Randomize.hh"
//#include "G4MTHepRandom.hh"

#include "CLHEP/Units/SystemOfUnits.h" // temp

#include <vector>
using std::vector;

class Level
{
public:
    Level(G4double, const vector<G4double>&, const vector<G4double>&, G4bool);
    ~Level();

    G4double Decay();

    G4double GetEnergy() const
    {
        return m_energy;
    }

    G4int GetNumberOfDaughters() const
    {
        return m_numberOfDaughters;
    }

    G4double GetDaughterEnergy(G4int i) const;

    G4bool IsEndState() const
    {
        return m_endState;
    }

private:
    G4double m_energy;
    G4int m_numberOfDaughters;
    G4bool m_endState;

    G4RandGeneral* ran = nullptr;

    G4double* m_daughterEnergies = nullptr;
    G4double* m_decayProbabilities = nullptr;
};

#endif // Level_hh
