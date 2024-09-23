// Levels.cc

#include "generator/GammaDecayScheme/Level.hh"


Level::Level(G4double energy, const vector<G4double> &daughterEnergies, const vector<G4double> &decayProbabilities, G4bool endState)
    : m_energy(energy), m_numberOfDaughters(0), m_endState(endState)
{
    if (decayProbabilities.size() != daughterEnergies.size())
    {
        G4cerr << "Error: number of daughter energies does not match number of decay probabilites" << G4endl;
        G4cerr << "Level " << energy << G4endl;
        G4cerr << daughterEnergies.size() << " daughters." << G4endl;
        G4cerr << decayProbabilities.size() << " decay probabilities." << G4endl;
        exit(99);
    }

    m_numberOfDaughters = daughterEnergies.size();

    if(!m_endState)
    {

        if (m_numberOfDaughters == 0)
        {
            G4cerr << "Error: no daughters given for non-end state." << G4endl;
            exit(99);
        }
        m_daughterEnergies = new G4double[m_numberOfDaughters];
        m_decayProbabilities = new G4double[m_numberOfDaughters];

        for(G4int i=0; i < m_numberOfDaughters; i++)
        {
            if (daughterEnergies.at(i) >= m_energy)
            {
                G4cerr << "Error: found daughter with higher energy than original state!" << G4endl;
                exit(99);
            }
            m_daughterEnergies[i] = daughterEnergies.at(i);
            m_decayProbabilities[i] = decayProbabilities.at(i);
        }
        ran = new G4RandGeneral(m_decayProbabilities, m_numberOfDaughters, 1);
    }
    else
    {
        if (m_numberOfDaughters > 0)
        {
            G4cerr << "Error: daughters given for end state." << G4endl;
            exit(99);
        }
    }
}

Level::~Level()
{
    //G4cout << "~Level" << G4endl;
    delete[] m_daughterEnergies;
    delete[] m_decayProbabilities;

    delete ran;
}

G4double Level::GetDaughterEnergy(G4int i) const
{
    if (i >= m_numberOfDaughters)
    {
        G4cerr << "Daughter energy out of range." << G4endl;
    }
    return m_daughterEnergies[i];
}

G4double Level::Decay()
{
    if(m_endState)
    {
        G4cerr << "Requested an end state to decay." << G4endl;
        exit(99);
    }

    return GetDaughterEnergy(G4int( ran->shoot() * m_numberOfDaughters + 0.5));
}
