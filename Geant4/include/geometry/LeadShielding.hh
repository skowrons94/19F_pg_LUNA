#ifndef LeadShielding_hh
#define LeadShielding_hh

#include "geometry/GeometryObject.hh"

class G4UIcmdWithAString;

#include <memory>
using std::shared_ptr;

class LeadShielding : public GeometryObject
{
public:
    LeadShielding();
    virtual ~LeadShielding() {}

    G4VPhysicalVolume *Construct();
    void ConstructSDandField() {}

    virtual void SetNewValue(G4UIcommand* command, G4String value);

private:
    shared_ptr<G4UIcmdWithAString> fSetupCmd;
    shared_ptr<G4UIcmdWithADoubleAndUnit> f55DegreePositionCmd;

    bool fSetupChosen = false;
    enum {setupBGO, setupGermanium0deg, setupGermanium55deg} fSetup;

    double f55degreePosition = 0.0;
};

#endif
