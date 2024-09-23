#ifndef TargetHolder_hh
#define TargetHolder_hh

#include "geometry/GeometryObject.hh"
class G4UIcmdWithAString;

#include <memory>
using std::shared_ptr;

class TargetHolder : public GeometryObject {

    public:
        TargetHolder();
        virtual ~TargetHolder();

        G4VPhysicalVolume *Construct();
        void ConstructSDandField() {};

        virtual void SetNewValue(G4UIcommand* command, G4String value);

    private:
        shared_ptr<G4UIcmdWithAString> fTargetCmd;

        bool fTargetChosen = false;
        enum {targetGraphite, targetSource, targetEvaporated} fTarget;
};

#endif
