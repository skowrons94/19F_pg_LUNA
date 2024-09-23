#include "GeometryManagerMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "GeometryManager.hh"

#include <memory>
using std::make_shared;

GeometryManagerMessenger::GeometryManagerMessenger(GeometryManager *geometryManager) :
    m_geometryManager(geometryManager)
{
    m_geometryDir = make_shared<G4UIdirectory>("/Geometry/");
    m_geometryDir->SetGuidance("Geometry control");

    m_addGeometryCmd = make_shared<G4UIcmdWithAString>("/Geometry/add", this);
    m_addGeometryCmd->SetGuidance("Enable geometry by name.");
    m_addGeometryCmd->SetParameterName("name", false);
    m_addGeometryCmd->SetToBeBroadcasted(false);

    m_ignoreOverlapsCmd = make_shared<G4UIcmdWithABool>("/Geometry/ignoreOverlaps", this);
    m_ignoreOverlapsCmd->SetGuidance("Do not stop when overlaps are found.");
    m_ignoreOverlapsCmd->SetParameterName("ignoreOverlaps", false);
    m_ignoreOverlapsCmd->SetToBeBroadcasted(false);

    m_doNotCheckOverlapsCmd = make_shared<G4UIcmdWithABool>("/Geometry/doNotCheckOverlaps", this);
    m_doNotCheckOverlapsCmd->SetGuidance("Do not check for overlaps .");
    m_doNotCheckOverlapsCmd->SetParameterName("ignoreOverlaps", false);
    m_doNotCheckOverlapsCmd->SetToBeBroadcasted(false);
}

void GeometryManagerMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if(command == m_addGeometryCmd.get())
    {
        m_geometryManager->AddGeometry(newValue);
    }
    else if (command == m_ignoreOverlapsCmd.get())
    {
        m_geometryManager->SetIgnoreOverlaps(m_ignoreOverlapsCmd->GetNewBoolValue(newValue));
    }
    else if (command == m_doNotCheckOverlapsCmd.get())
    {
        m_geometryManager->SetDoNotCheckOverlaps(m_doNotCheckOverlapsCmd->GetNewBoolValue(newValue));
    }
    else
    {
        G4cerr << "Unknown command in GeometryManagerMessenger!" << G4endl;
    }
}
