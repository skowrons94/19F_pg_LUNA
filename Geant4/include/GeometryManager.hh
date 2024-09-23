#ifndef GeometryManager_hh
#define GeometryManager_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4String.hh"
#include "geometry/GeometryObject.hh"

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <memory>
using std::shared_ptr;

class GeometryManagerMessenger;
class GeometryObject;

/// Geometry Manager Class

/** This is a singleton class that manages the GeometryObjects. All objects are
  * instantiated in the constructor, but only those enabled by macro commands
  * are being constructed when the run is initialized.
  *
  * Macro commands are defined in GeometryManagerMessenger class. See
  * documentation there on how to enable geometries and change options for
  * checking overlaps.
  *
**/

class GeometryManager : public G4VUserDetectorConstruction
{
private:
    // Private constructor for singleton pattern
    GeometryManager();
    static GeometryManager* instance;

public:
    ~GeometryManager() {}

    static GeometryManager* GetInstance()
    {
        if (GeometryManager::instance == nullptr)
            new GeometryManager();
        return GeometryManager::instance;
    }

    /// Select a geometry object by name to be built
    void AddGeometry(const G4String &name);

    /// Check whether a geometry object is in the list to be built
    G4bool HasGeometry(const G4String &name);

    /// Get pointer to geometry object by name
    GeometryObject *GetGeometry(const G4String &name);

    /// Get value of option to ignore (not throw an exception) in case of overlaps
    G4bool GetIgnoreOverlaps()
    {
        return m_ignoreOverlaps;
    }

    /// Set value of option to ignore (not throw an exception) in case of overlaps
    void SetIgnoreOverlaps(G4bool newValue)
    {
        m_ignoreOverlaps = newValue;
    }

    /// Get value of option to not check for overlaps
    G4bool GetDoNotCheckOverlaps()
    {
        return m_doNotCheckOverlaps;
    }

    /// Set value of option to not check for overlaps
    void SetDoNotCheckOverlaps(G4bool newValue)
    {
        m_doNotCheckOverlaps = newValue;
    }

private:
    /// Define all materials
    void ConstructMaterials();

public:
    /// Construct geometry
    G4VPhysicalVolume *Construct();

    /// Sensitive detectors and fields
    void ConstructSDandField();

    /// Set up outputs
    void SetupOutputs();

    /// Let geometry objects fill their outputs
    void FillOutputs(const G4Event *event);

private:

    shared_ptr<GeometryManagerMessenger> m_GeometryManagerMessenger;

    vector<shared_ptr<GeometryObject>> m_allGeometries, m_selectedGeometries;

    G4bool m_ignoreOverlaps = true;
    G4bool m_doNotCheckOverlaps = false;
};

#endif // GeometryManager_hh
