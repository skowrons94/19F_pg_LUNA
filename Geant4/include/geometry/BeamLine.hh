#ifndef BeamLine_hh
#define BeamLine_hh

#include "geometry/GeometryObject.hh"

/**
    Approximate BeamLine Geometry

    Mostly used for neutron capture background simulations. Dimensions are
    "educated guesses". If the precise geometry of this part is expected to be
    of importance for the simulation results, the dimensions should be checked
    carefully!
**/

class BeamLine : public GeometryObject
{
public:
    BeamLine() : GeometryObject("BeamLine") {}
    virtual ~BeamLine() {}

    G4VPhysicalVolume *Construct();
    void ConstructSDandField() {};
};

#endif
