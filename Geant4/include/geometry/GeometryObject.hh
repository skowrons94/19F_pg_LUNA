#ifndef GeometryObject_hh
#define GeometryObject_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4UImessenger.hh"
#include "G4VisAttributes.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

#include "GeometryManager.hh"

#include <memory>
using std::shared_ptr;

class G4PVPlacement;

/// Base class for Geometry Objects

/** Geometry objects, the building blocks of the experimental geometry in this
  * framework derive from this class.
  *
  * Each geometry object has a mother volume, in which it will be placed.
  *
  * The position and rotation of the object can be set by automatically defined
  * commands:
  *
  *   /Geometry/objectName/position x y z unit
  *
  *   /Geometry/objectName/rotateX angle deg
  *   /Geometry/objectName/rotateY angle deg
  *   /Geometry/objectName/rotateZ angle deg
  *
  * Especially in view of rotations, the point of origin for each object should
  * be chosen carefully. For a HPGe or other "can style" detector, for example,
  * the center of the front face can be a good choice for convenient rotation
  * and placement.
  * Target chamber type objects may benefit from having the center of the target
  * surface at (0,0,0).
  *
  * Apart from the geometry itself, a GeometryObject also provides virtual
  * methods to define and to fill outputs of the analysis class.
  *
  * The SetupOutput() has to be called from the BeginOfRunAction.
  *
  * FillOutput() is called at the end of each event.
  *
**/

class GeometryObject : public G4VUserDetectorConstruction, public G4UImessenger {

    public:
        GeometryObject(const G4String &name);
        virtual ~GeometryObject() {}

        const G4String &GetName() const {return m_name;}
        G4bool GetIsSensitive() const {return m_buildSensitiveVolumes;}

        /// Purely virtual functions, to be implemented by child classes
        virtual G4VPhysicalVolume *Construct() = 0;
        virtual void ConstructSDandField() = 0;

        void SetMotherVolume(G4LogicalVolume *volume)
        {
            m_motherVolume = volume;
        }

        G4LogicalVolume* GetMotherVolume() const {return m_motherVolume;}

        const G4ThreeVector &GetPosition() const {return m_position;}
        G4RotationMatrix *GetRotation() {return m_rotation;}

        G4Transform3D GetTransform3D(G4RotationMatrix ownRotation, G4ThreeVector relativePosition);

        G4PVPlacement *PlaceVolume(G4LogicalVolume *logicalVolume, G4LogicalVolume *motherVolume = nullptr, G4ThreeVector position = G4ThreeVector(), G4RotationMatrix rotation = G4RotationMatrix(), int copyNr = 0, bool noExternalRotation = false);
        G4PVPlacement *PlaceVolumeInternal(G4LogicalVolume *logicalVolume, G4LogicalVolume *motherVolume = nullptr, G4ThreeVector position = G4ThreeVector(), G4RotationMatrix rotation = G4RotationMatrix(), int copyNr = 0);

        virtual void SetupOutput() {}
        virtual void FillOutput(const G4Event */*event*/) {}
        virtual void FillPosition(G4ThreeVector &) {};
  
        virtual void SetNewValue(G4UIcommand* command, G4String value);

    private:
        G4String m_name;

        G4ThreeVector m_position;
        G4RotationMatrix *m_rotation;
        G4bool m_buildSensitiveVolumes = true;

        G4LogicalVolume *m_motherVolume;

        G4String m_cmdDirName;

        shared_ptr<G4UIdirectory> m_cmdDir;

        shared_ptr<G4UIcmdWith3VectorAndUnit> m_cmdSetPosition;
        shared_ptr<G4UIcmdWithADoubleAndUnit> m_cmdRotateX;
        shared_ptr<G4UIcmdWithADoubleAndUnit> m_cmdRotateY;
        shared_ptr<G4UIcmdWithADoubleAndUnit> m_cmdRotateZ;

        shared_ptr<G4UIcmdWithABool> m_cmdSensitiveVolumes;
};

#endif
