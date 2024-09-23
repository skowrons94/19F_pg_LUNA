#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "G4UImessenger.hh"
class G4UIcmdWithAString;

#include <memory>
using std::shared_ptr;


/// Run action class
///

class RunAction : public G4UserRunAction, public G4UImessenger
{
public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    virtual void SetNewValue(G4UIcommand*, G4String);

private:
    G4String fFileName = "";

    shared_ptr<G4UIcmdWithAString> m_setFileNameCmd;
};


#endif

