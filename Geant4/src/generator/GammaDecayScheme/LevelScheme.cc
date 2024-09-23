#include "generator/GammaDecayScheme/LevelScheme.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using CLHEP::keV;

#include <fstream>
#include <sstream>
#include <vector>

using std::ifstream;
using std::stringstream;
using std::vector;

#include <memory>
using std::make_shared;

#include <string>
using std::string;

#include <stdexcept>
using std::runtime_error;

LevelScheme::LevelScheme(const G4String &levelFileName)
{
    ifstream file;
    file.open(levelFileName.c_str());

    if (!file)
    {
        G4cerr << "Could not open level file \"" << levelFileName << "\"!" << G4endl;
        throw runtime_error("Failed to open level file in LevelScheme::LevelScheme().");
    };

    G4double highestState = -1 * keV;
    vector<G4double> p, d;

    G4double energy = -1*keV;
    G4String line, word;
    stringstream ss;

    enum ReadStatus {reset, idle, state, daughters, probabilities, finishedlevel};
    ReadStatus status = idle;

    G4bool endState = false;

    //read our file three lines at a time
    while(getline(file, line))
    {

        if (status == reset)
        {
            energy = -1*keV;
            endState = false;
            p.clear();
            d.clear();
            status = idle;
        }

        ss.clear();
        ss.str(line);

        // is this line empty?
        if (!(ss >> word))
            continue;

        // is this line a (pure) comment?
        if (word.find("//") == 0)
            continue;

        // line should be interesting, let's have a look

        ss.clear();
        ss.str(line);

        switch (status)
        {
        case idle:
            ss >> word;
            if (word != "State")
            {
                G4cerr << "Expected to read 'State', found '" << word << "'." << G4endl;
                exit(99);
            }
            ss >> word;
            energy = atof(word.c_str()) * keV;
            status = daughters;
            break;

        case daughters:
            ss >> word;
            if (word == "end")
            {
                endState = true;
                status = finishedlevel;
                break;
            }
            ss.clear();
            ss.str(line);
            while (ss >> word)
            {
                if (line.find("//") == 0)
                    break;
                d.push_back(atof(word.c_str()) * keV);
            }
            status = probabilities;
            break;

        case probabilities:
            while (ss >> word)
            {
                if (line.find("//") == 0)
                    break;
                p.push_back(atof(word.c_str()));
            }
            status = finishedlevel;
            break;

        case reset:
        case finishedlevel:
        default:
            G4cerr << "Unexpected state while reading the level file: " << state << G4endl;
            exit(99);
            break;
        }

        if (status == finishedlevel)
        {

            if (energy > highestState)
                highestState = energy;

            NewLevel(energy, d, p, endState);

            status = reset;
        }
    }

    file.close();

    SetStartLevel(highestState);
}


void LevelScheme::NewLevel(G4double energy, const vector<G4double> &daughterEnergies, const vector<G4double> &decayProbabilities, G4bool endState)
{
    if (levels.find(energy) != levels.end())
    {
        G4cerr << "Error trying to insert level at " << energy / keV << " keV. Level exists already." << G4endl;
        throw runtime_error("LevelScheme::NewLevel(): Tried to insert level at existing energy.");
    }
    //G4cout << "Inserted level at " << energy / keV << "keV" << G4endl;
    levels[energy] = make_shared<Level>(energy, daughterEnergies, decayProbabilities, endState);
}

void LevelScheme::SetLevel(G4double energy)
{
    const auto it = levels.find(energy);
    if (it == levels.end())
    {
        G4cerr << "Error trying to find level at " << energy / keV << " keV. Level doesn't exist." << G4endl;
        throw runtime_error("LevelScheme::SetLevel(): Could not find level.");
    }
    currentLevel = it->second;
}

void LevelScheme::SetStartLevel(G4double energy)
{
    const auto it = levels.find(energy);
    if (it == levels.end())
    {
        G4cerr << "Error trying to find level at " << energy / keV << " keV. Level doesn't exist." << G4endl;
        throw runtime_error("LevelScheme::SetStartLevel(): Could not find level.");
    }
    startLevel = it->second;
    CheckConsistency(energy);
}

void LevelScheme::Start()
{
    currentLevel = startLevel;
}

G4double LevelScheme::Decay()
{
    G4double currentEnergy = currentLevel->GetEnergy();
    G4double daughterEnergy = currentLevel->Decay();
    SetLevel(daughterEnergy);

    return currentEnergy - daughterEnergy;
}


void LevelScheme::CheckConsistency(G4double startEnergy)
{
    const auto it = levels.find(startEnergy);
    if (it == levels.end())
    {
        G4cerr << "Consistency Check failed: level " << startEnergy / keV << " not found." << G4endl;
        throw runtime_error("LevelScheme::CheckConsistency(): Consistency check failed.");
    }
    auto level = it->second;

    if (level->IsEndState())
    {
        return;
    }

    auto n = level->GetNumberOfDaughters();
    for (G4int i = 0; i < n; i++)
    {
        CheckConsistency(level->GetDaughterEnergy(i));
    }
}
