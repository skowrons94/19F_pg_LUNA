#ifndef LevelScheme_hh
#define LevelScheme_hh

#include "globals.hh"

#include <map>
using std::map;

#include <vector>
using std::vector;

#include <memory>
using std::shared_ptr;

#include "generator/GammaDecayScheme/Level.hh"


class LevelScheme {
    public:
        LevelScheme(const G4String &levelFileName);
        ~LevelScheme() {}

        void NewLevel(G4double energy, const vector<G4double> &daughterEnergies, const vector<G4double> &decayProbabilities, G4bool endState);
        void SetLevel(G4double energy);
        void SetStartLevel(G4double energy);
        void Start();

        G4bool IsAtEndState() const {return currentLevel->IsEndState();}
        G4double Decay();

    private:
        void CheckConsistency(G4double startEnergy);
        map<G4double, shared_ptr<Level>> levels;
        shared_ptr<Level> startLevel, currentLevel;

};

#endif // LevelScheme_hh
