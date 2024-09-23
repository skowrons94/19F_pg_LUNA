#include "GeometryManager.hh"
#include "GeometryManagerMessenger.hh"

#include "geometry/GeometryObject.hh"

#include "geometry/LeadShielding.hh"
#include "geometry/BGODetector.hh"
#include "geometry/TargetChamberMg25.hh"
#include "geometry/TargetChamberO17.hh"
#include "geometry/TargetHolder.hh"
#include "geometry/BeamLine.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4Box.hh"
#include "G4Event.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <stdexcept>
using std::runtime_error;

#include <memory>
using std::make_shared;

using CLHEP::g;
using CLHEP::mm;
using CLHEP::cm;

using CLHEP::cm3;
using CLHEP::perCent;
using CLHEP::deg;

GeometryManager* GeometryManager::instance = nullptr;

GeometryManager::GeometryManager() :
    m_GeometryManagerMessenger(make_shared<GeometryManagerMessenger>(this))
{
    if (GeometryManager::instance != nullptr)
    {
        G4cout << "GeometryManager() called more than once. This should not happen." << G4endl;
        exit(99);
    }

    /// Fill the list of possible geometry objects
    m_allGeometries.push_back(make_shared<BGODetector>());
    m_allGeometries.push_back(make_shared<LeadShielding>());
    m_allGeometries.push_back(make_shared<TargetChamberMg25>());
    m_allGeometries.push_back(make_shared<TargetChamberO17>());
    m_allGeometries.push_back(make_shared<TargetHolder>());
    m_allGeometries.push_back(make_shared<BeamLine>());
    // ...

    GeometryManager::instance = this;
}


void GeometryManager::AddGeometry(const G4String &name)
{
    // check wether this geometry is already in the list of objects to build
    for (const auto &geo : m_selectedGeometries)
    {
        if (geo->GetName() == name)
        {
            G4cout << "Tried to add geometry " << name << ", but was already present." << G4endl;
            exit(99);
        }
    }

    // try to add this geometry
    for (const auto &geo : m_allGeometries)
    {
        if (geo->GetName() == name)
        {
            m_selectedGeometries.push_back(geo);
            return;
        }
    }

    G4cout << "Tried to add geometry " << name << ", but wasn't found." << G4endl;
    exit(99);
}

G4bool GeometryManager::HasGeometry(const G4String &name)
{
    return !(GetGeometry(name) == nullptr);
}


GeometryObject *GeometryManager::GetGeometry(const G4String &name)
{
    // check wether this geometry is already in the list of objects to build
    for (const auto &geo : m_selectedGeometries)
    {
        if (geo->GetName() == name)
            return geo.get();
    }
    return nullptr;
}

void GeometryManager::ConstructMaterials()
{
    G4NistManager* man = G4NistManager::Instance();

    // Collect pure materials
    G4Material* matFe = man->FindOrBuildMaterial("G4_Fe");
    G4Material* matCr = man->FindOrBuildMaterial("G4_Cr");
    G4Material* matNi = man->FindOrBuildMaterial("G4_Ni");
    G4Material* matMo = man->FindOrBuildMaterial("G4_Mo");

    G4Material* matCu = man->FindOrBuildMaterial("G4_Cu");
    G4Material* matZn = man->FindOrBuildMaterial("G4_Zn");
    G4Material* matPb = man->FindOrBuildMaterial("G4_Pb");
    
    G4Material* matMg = man->FindOrBuildMaterial("G4_Mg");
    G4Material* matMn = man->FindOrBuildMaterial("G4_Mn");
    G4Material* matSi = man->FindOrBuildMaterial("G4_Si");
    G4Material* matTi = man->FindOrBuildMaterial("G4_Ti");
    G4Material* matAl = man->FindOrBuildMaterial("G4_Al");

    // Steel-316: Stainless Steel
    G4Material* matSteel = new G4Material("Steel-316", 8.0*g/cm3, 4);
    matSteel->AddMaterial(matFe, 69*perCent);
    matSteel->AddMaterial(matCr, 17*perCent);
    matSteel->AddMaterial(matNi, 12*perCent);
    matSteel->AddMaterial(matMo,  2*perCent);

    G4Material* matBrass = new G4Material("Brass", 8.52*g/cm3, 3);
    matBrass->AddMaterial(matCu, 65*perCent);
    matBrass->AddMaterial(matZn, 33*perCent);
    matBrass->AddMaterial(matPb, 2*perCent);

    G4Material* matAlloy = new G4Material("Alloy", 2.70*g/cm3, 9);
    matAlloy->AddMaterial(matAl, 98.0*perCent);
    matAlloy->AddMaterial(matCr, 0.05*perCent);
    matAlloy->AddMaterial(matCu, 0.10*perCent);
    matAlloy->AddMaterial(matFe, 0.30*perCent);
    matAlloy->AddMaterial(matMg, 0.60*perCent);
    matAlloy->AddMaterial(matMn, 0.10*perCent);
    matAlloy->AddMaterial(matSi, 0.60*perCent);
    matAlloy->AddMaterial(matZn, 0.15*perCent);
    matAlloy->AddMaterial(matTi, 0.10*perCent);
    
}

G4VPhysicalVolume *GeometryManager::Construct()
{

    G4cout << "GeometryManager::Construct()" << G4endl;

    //
    ConstructMaterials();

    // construct world volume
    G4NistManager* man = G4NistManager::Instance();
    G4Material* matAir = man->FindOrBuildMaterial("G4_AIR");
    G4Material* matVac = man->FindOrBuildMaterial("G4_Galactic");

    //create world and experiment hall

    G4double worldSize = 5*CLHEP::m;

    G4Box* worldSolid = new G4Box("worldSolid", 0.5*worldSize, 0.5*worldSize, 0.5*worldSize);
    G4LogicalVolume* worldLog = new G4LogicalVolume(worldSolid, matVac, "worldLog", 0, 0, 0);
    //worldLog->SetVisAttributes(G4VisAttributes::Invisible);

    // construct geometry objects
    for (auto &geo : m_selectedGeometries)
    {
        G4cout << "calling " << geo->GetName() << "::Construct()" << G4endl;
        geo->SetMotherVolume(worldLog);
        geo->Construct();
    }

    G4VPhysicalVolume* worldPhys = new G4PVPlacement(nullptr, G4ThreeVector(), worldLog, "World_phys",
            0, false, 0);

    return worldPhys;
}

void GeometryManager::ConstructSDandField()
{
    for (auto &geo : m_selectedGeometries)
    {
        if (geo->GetIsSensitive())
        {
            geo->ConstructSDandField();
        }
    }
}

void GeometryManager::SetupOutputs()
{
    for (auto &geo : m_selectedGeometries)
    {
        if (geo->GetIsSensitive())
        {
            geo->SetupOutput();
        }
    }
}

void GeometryManager::FillOutputs(const G4Event *event)
{
    for (auto &geo : m_selectedGeometries)
    {
        if (geo->GetIsSensitive())
        {
            geo->FillOutput(event);
        }
    }
}

