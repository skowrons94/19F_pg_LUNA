#ifndef TargetChamberMg25_hh
#define TargetChamberMg25_hh

#include "geometry/GeometryObject.hh"

class TargetChamberMg25 : public GeometryObject {

    public:
        TargetChamberMg25();
        virtual ~TargetChamberMg25();

        G4VPhysicalVolume *Construct();
        void ConstructSDandField() {};

        virtual void SetNewValue(G4UIcommand* command, G4String value);

    private:
};

#endif
