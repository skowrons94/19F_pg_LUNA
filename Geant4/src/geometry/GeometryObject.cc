#include "geometry/GeometryObject.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"

#include <memory>
using std::make_shared;

GeometryObject::GeometryObject(const G4String &name) :
    m_name(name),
    m_position(0, 0, 0),
    m_rotation(new G4RotationMatrix()),
    m_motherVolume(NULL),
    m_cmdDirName("/Geometry/" + name + "/")
{
    m_cmdDir = make_shared<G4UIdirectory>(m_cmdDirName);
    m_cmdDir->SetGuidance(("commands for " + m_name + " geometry.").c_str());

    m_cmdSetPosition = make_shared<G4UIcmdWith3VectorAndUnit>((m_cmdDirName + "position").c_str(), static_cast<G4UImessenger*>(this));
    m_cmdSetPosition->SetGuidance(("Position of " + m_name + " .").c_str());
    m_cmdSetPosition->SetParameterName("X", "Y", "Z", false);
    m_cmdSetPosition->SetUnitCategory("Length");

    m_cmdRotateX = make_shared<G4UIcmdWithADoubleAndUnit>((m_cmdDirName + "rotateX").c_str(), static_cast<G4UImessenger*>(this));
    m_cmdRotateX->SetGuidance(("Rotate " + m_name + " around x-axis.").c_str());
    m_cmdRotateX->SetParameterName("angĺeX", false);
    m_cmdRotateX->SetUnitCategory("Angle");

    m_cmdRotateY = make_shared<G4UIcmdWithADoubleAndUnit>((m_cmdDirName + "rotateY").c_str(), static_cast<G4UImessenger*>(this));
    m_cmdRotateY->SetGuidance(("Rotate " + m_name + " around y-axis.").c_str());
    m_cmdRotateY->SetParameterName("angĺeY", false);
    m_cmdRotateY->SetUnitCategory("Angle");

    m_cmdRotateZ = make_shared<G4UIcmdWithADoubleAndUnit>((m_cmdDirName + "rotateZ").c_str(), static_cast<G4UImessenger*>(this));
    m_cmdRotateZ->SetGuidance(("Rotate " + m_name + " around z-axis.").c_str());
    m_cmdRotateZ->SetParameterName("angĺeZ", false);
    m_cmdRotateZ->SetUnitCategory("Angle");

    m_cmdSensitiveVolumes = make_shared<G4UIcmdWithABool>((m_cmdDirName + "sensitiveVolumes").c_str(), static_cast<G4UImessenger*>(this));
    m_cmdSensitiveVolumes->SetGuidance("Set sensitive detector volumes enabled/disabled (= true/false).");
}

G4Transform3D GeometryObject::GetTransform3D(G4RotationMatrix ownRotation, G4ThreeVector relativePosition)
{
    return G4Transform3D(((*GetRotation())*ownRotation), GetPosition() + (*GetRotation())(relativePosition));
}

G4PVPlacement* GeometryObject::PlaceVolume(G4LogicalVolume *logicalVolume, G4LogicalVolume *motherVolume, G4ThreeVector position, G4RotationMatrix rotation, int copyNr, bool noExternalRotation)
{
    if (motherVolume == nullptr)
    {
        motherVolume = GetMotherVolume();
    }

    G4PVPlacement *phyVol = nullptr;

    /// FIXME: Something's odd here, these two cases shouldn't be identical?!
    if (noExternalRotation)
    {
        phyVol =
            new G4PVPlacement(G4Transform3D(rotation, position),
                              logicalVolume, logicalVolume->GetName(),
                              motherVolume, false, copyNr, false);
    }
    else
    {
        phyVol =
            new G4PVPlacement(GetTransform3D(rotation, position),
                              logicalVolume, logicalVolume->GetName(),
                              motherVolume, false, copyNr, false);
    }

    // check for overlaps
    if (!GeometryManager::GetInstance()->GetDoNotCheckOverlaps())
    {
        if (phyVol->CheckOverlaps())
        {
            G4cerr << "Overlaps found!" << G4endl;
            if (!GeometryManager::GetInstance()->GetIgnoreOverlaps())
                exit(99);
        }
    }
    return phyVol;
}

G4PVPlacement* GeometryObject::PlaceVolumeInternal(G4LogicalVolume *logicalVolume, G4LogicalVolume *motherVolume, G4ThreeVector position, G4RotationMatrix rotation, int copyNr)
{
    return PlaceVolume(logicalVolume, motherVolume, position, rotation, copyNr, true);
}


void GeometryObject::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == m_cmdSetPosition.get())
    {
        m_position = static_cast<G4UIcmdWith3VectorAndUnit*>(command)->GetNew3VectorValue(newValue);
    }
    else if (command == m_cmdRotateX.get())
    {
        m_rotation->rotateX(static_cast<G4UIcmdWithADoubleAndUnit*>(command)->GetNewDoubleValue(newValue));
    }
    else if (command == m_cmdRotateY.get())
    {
        m_rotation->rotateY(static_cast<G4UIcmdWithADoubleAndUnit*>(command)->GetNewDoubleValue(newValue));
    }
    else if (command == m_cmdRotateZ.get())
    {
        m_rotation->rotateZ(static_cast<G4UIcmdWithADoubleAndUnit*>(command)->GetNewDoubleValue(newValue));
    }
    else if (command == m_cmdSensitiveVolumes.get())
    {
        m_buildSensitiveVolumes = static_cast<G4UIcmdWithABool*>(command)->GetNewBoolValue(newValue);
    }
    else
    {
        G4cerr << "Unknown command in GeometryObject!" << G4endl;
    }
}
