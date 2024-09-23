#include "RunAction.hh"
#include "Analysis.hh"
#include "GeometryManager.hh"

#include "G4UIcmdWithAString.hh"

#include <memory>
using std::make_shared;

RunAction::RunAction() : G4UserRunAction()
{
    m_setFileNameCmd = make_shared<G4UIcmdWithAString>("/run/setFileName", this);
    m_setFileNameCmd->SetGuidance("Set output file name.");
    m_setFileNameCmd->SetParameterName("file name", false);
}

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    // Inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

    // Set analysis manager options
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetNtupleMerging(true);

    // Open an output file
    //
    analysisManager->OpenFile(fFileName);

    // Setup outputs of the geometry objects
    auto *gm = GeometryManager::GetInstance();
    gm->SetupOutputs();
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    // Save and close output
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

void RunAction::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == m_setFileNameCmd.get())
    {
        fFileName = newValue;
    }
    else
    {
        G4cerr << "Unknown command in RunAction!" << G4endl;
    }
}
