/** \brief LUNA BL2 Detector Simulation

    \details A collection of geometries and useful tools/classes to simulate
             the solid target station on beam line 2 at LUNA-400.
**/

#include "GeometryManager.hh"

#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "PhysicsList.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


void PrintUsage()
{
    G4cerr << " Usage: " << G4endl;
    G4cerr << " G4BL2 [-m macro ] [-u UIsession]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    /// Parse command line arguments
    //

    G4String macro = "";
    G4String session = "";

    {
        if ( argc > 5 )
        {
            PrintUsage();
            return 1;
        }

        for (G4int i = 1; i < argc; i = i + 2)
        {
            if (G4String(argv[i]) == "-m")
            {
                macro = argv[i+1];
            }
            else if (G4String(argv[i]) == "-u")
            {
                session = argv[i+1];
            }
            else
            {
                PrintUsage();
                return 1;
            }
        }
    }

    /// Interactive mode if no macro is provided
    //
    G4UIExecutive* ui = nullptr;
    if (!macro.size())
    {
        ui = new G4UIExecutive(argc, argv, session);
    }

    /// Choose the Random engine
    //
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

    /// Construct the run manager
    //
#ifdef G4MULTITHREADED
    auto runManager = new G4MTRunManager();
#else
    auto runManager = new G4RunManager;
#endif

    /// Set mandatory initialization classes
    //
    {
        // Implementation of G4VUserDetectorConstruction
        runManager->SetUserInitialization(GeometryManager::GetInstance());

        // Physics list
	//	G4PhysListFactory factory;
	//	G4VModularPhysicsList* physList = factory.GetReferencePhysList("FTFP_BERT_LIV");
	//	runManager->SetUserInitialization(physList);
	runManager->SetUserInitialization(new PhysicsList);

        // Implementation of G4VUserActionInitialization
        auto actionInitialization = new ActionInitialization();
        runManager->SetUserInitialization(actionInitialization);
    }

    /// Initialize visualization
    //
    auto visManager = new G4VisExecutive();
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    /// Get the pointer to the User Interface manager
    //
    auto UImanager = G4UImanager::GetUIpointer();

    /// Process macro or start UI session
    //
    if (macro.size())
    {
        // batch mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macro);
    }
    else
    {
        // interactive mode : define UI session
        UImanager->ApplyCommand("/control/execute mac/init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

    //delete visManager; // uncommenting this causes seg-fault when closing GUI in interactive mode(?)
    delete runManager;

    return 0;
}
