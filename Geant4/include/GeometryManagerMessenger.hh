#ifndef GeometryManagerMessenger_hh
#define GeometryManagerMessenger_hh

#include "G4UImessenger.hh"

#include <memory>
using std::shared_ptr;

class G4String;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class GeometryManager;

/// Messenger Class for the GeometryManager

/** Defines the following commands:
  *
  *   /Geometry/add ObjectName
  *      ... adds an element to the list of objects to be built
  *
  *   /Geometry/ignoreOverlaps True/False
  *      ... if True, do not throw exception when finding overlaps
  *
  *   /Geometry/doNotCheckOverlaps True/False
  *      ... if True, do not check for overlapping geometries
**/

class GeometryManagerMessenger : public G4UImessenger {

    public:
        GeometryManagerMessenger(GeometryManager*);
        virtual ~GeometryManagerMessenger() {}

        void SetNewValue(G4UIcommand*, G4String);

    private:
        GeometryManager* m_geometryManager;

        shared_ptr<G4UIdirectory>      m_geometryDir;

        shared_ptr<G4UIcmdWithAString> m_addGeometryCmd;

        shared_ptr<G4UIcmdWithABool>   m_ignoreOverlapsCmd,
                                       m_doNotCheckOverlapsCmd;
};




#endif // GeometryManagerMessenger_hh
