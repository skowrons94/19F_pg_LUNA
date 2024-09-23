#ifndef BGODetector_hh
#define BGODetector_hh

#include "geometry/GeometryObject.hh"

#include "G4ThreeVector.hh"

#include <memory>
using std::shared_ptr;

class G4MultiFunctionalDetector;
class G4PSEnergyDeposit;

/**
    BGO Detector Geometry

    Dimensions based on Scionix' drawing of the detector. Original
    implementation by A. Best. Slight changes (PMT cutouts) after that.

    A few dimensions had to be "estimated" from the drawings. Volumes without a
    description in the drawing are assumed empty.
**/

class BGODetector : public GeometryObject
{
public:
    BGODetector() : GeometryObject("BGODetector") {}
    ~BGODetector() {}

    G4VPhysicalVolume *Construct();
    void ConstructSDandField();

    virtual void SetupOutput();
    virtual void FillOutput(const G4Event *event);

    void FillPosition(G4ThreeVector &vec);

private:
    int fHCID = -1; /// Hit collection ID for energy deposition scorer
    int fTupleID = -1; /// tuple ID for energy deposition output
};

#endif // BGODetector_hh
