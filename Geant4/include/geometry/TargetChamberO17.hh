#ifndef TargetChamberO17_hh
#define TargetChamberO17_hh

#include "geometry/GeometryObject.hh"

class TargetChamberO17 : public GeometryObject {

    public:
        TargetChamberO17();
        virtual ~TargetChamberO17();

        G4VPhysicalVolume *Construct();
        void ConstructSDandField() {};

        virtual void SetNewValue(G4UIcommand* command, G4String value);

    private:
};

#endif
