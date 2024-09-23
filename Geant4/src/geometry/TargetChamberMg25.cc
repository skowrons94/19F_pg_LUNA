#include "geometry/TargetChamberMg25.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"

#include "G4Polycone.hh"
#include "G4Tubs.hh"

#include "G4SubtractionSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using CLHEP::cm;
using CLHEP::mm;
using CLHEP::um;
using CLHEP::pi;
using CLHEP::deg;
using CLHEP::perCent;
using CLHEP::cm3;
using CLHEP::g;


TargetChamberMg25::TargetChamberMg25() : GeometryObject("TargetChamberMg25")
{
}

TargetChamberMg25::~TargetChamberMg25()
{
}

G4VPhysicalVolume *TargetChamberMg25::Construct()
{

    // No water, for now
    const double beamLineLength = 275.0*mm;
    const double flangeWidth = 12.5*mm; // ??

    G4NistManager* man = G4NistManager::Instance();

    G4Material* matSteel  = G4Material::GetMaterial("Steel-316");
    auto matTargetChamber = matSteel;

    G4Material* matPVC = man->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
    G4Material* matTa = man->FindOrBuildMaterial("G4_Ta");
    G4Material* matVacuum = man->FindOrBuildMaterial("G4_Galactic");

    // ----------------- LUNA beam line ----------------------------------------------------------------
    G4Tubs* holder_SS1 = new G4Tubs(GetName() + "_holderSS1", 11.95*mm, 25.00*mm,   2.5*mm, 0.0, 360.*deg);
    G4Tubs* holder_SS2 = new G4Tubs(GetName() + "_holderSS2", 20.94*mm, 25.00*mm,   1.0*mm, 0.0, 360.*deg);
    G4Tubs* holder_SS3 = new G4Tubs(GetName() + "_holderSS3", 23.98*mm, 25.00*mm,   1.5*mm, 0.0, 360.*deg);
    G4Tubs* holder_SS4 = new G4Tubs(GetName() + "_holderSS4", 23.98*mm, 29.48*mm,   5.0*mm, 0.0, 360.*deg);
    G4Tubs* beamline_tub    = new G4Tubs(GetName() + "_beamline",    23.5*mm,25.*mm,        0.5*beamLineLength-10*um, 0.0, 360.*deg);
    G4Tubs* beamlinevac_tub = new G4Tubs(GetName() + "_beamlinevac",       0,23.5*mm-10*um, 0.5*beamLineLength, 0.0, 360.*deg);
    G4Tubs* flange_tub = new G4Tubs(GetName() + "_flange", 25*mm + 10*um, 50.*mm, 0.5*flangeWidth, 0.0, 360.*deg);
    G4Tubs* backing_tub = new G4Tubs(GetName() + "_backing", 0.*mm, 19.*mm, 0.2*mm, 0., 2.*pi);
    G4Tubs* cooling_tub = new G4Tubs(GetName() + "_cooling", 0.*mm, 29.48*mm, 2.5*mm, 0., 2.*pi);

    G4LogicalVolume* holderSS1_log =
        new G4LogicalVolume(holder_SS1, matTargetChamber,
                            holder_SS1->GetName());
    G4LogicalVolume* holderSS2_log =
        new G4LogicalVolume(holder_SS2, matTargetChamber,
                            holder_SS2->GetName());
    G4LogicalVolume* holderSS3_log =
        new G4LogicalVolume(holder_SS3, matTargetChamber,
                            holder_SS3->GetName());
    G4LogicalVolume* holderSS4_log =
        new G4LogicalVolume(holder_SS4, matTargetChamber,
                            holder_SS4->GetName());
    G4LogicalVolume* beamline_log =
        new G4LogicalVolume(beamline_tub, matTargetChamber,
                            beamline_tub->GetName());

    G4LogicalVolume* beamlinevac_log =
        new G4LogicalVolume(beamlinevac_tub, matVacuum,
                            beamlinevac_tub->GetName());

    G4LogicalVolume* flange_log =
        new G4LogicalVolume(flange_tub, matTargetChamber,
                            flange_tub->GetName());
    G4LogicalVolume* backing_log =
        new G4LogicalVolume(backing_tub, matTa,
                            backing_tub->GetName());
    G4LogicalVolume* cooling_log =
        new G4LogicalVolume(cooling_tub, matPVC,
                            cooling_tub->GetName());

    PlaceVolume(holderSS1_log, GetMotherVolume(), -G4ThreeVector(0.,0.,2.5*mm));
    PlaceVolume(holderSS2_log, GetMotherVolume(), -G4ThreeVector(0.,0.,-1.*mm));
    PlaceVolume(holderSS3_log, GetMotherVolume(), -G4ThreeVector(0.,0.,-3.5*mm));
    PlaceVolume(holderSS4_log, GetMotherVolume(), -G4ThreeVector(0.,0.,-10.*mm));

    PlaceVolume(beamline_log,      GetMotherVolume(), -G4ThreeVector(0.,0.,0.5*beamLineLength + 5.*mm));
    PlaceVolume(beamlinevac_log,   GetMotherVolume(), -G4ThreeVector(0.,0.,0.5*beamLineLength + 5.*mm));

    PlaceVolume(flange_log,   GetMotherVolume(),   -G4ThreeVector(0.,0.,0.5*flangeWidth - 12.5*mm + 280*mm));
    PlaceVolume(cooling_log,    GetMotherVolume(), -G4ThreeVector(0.,0.,-17.505*mm));
    PlaceVolume(backing_log, GetMotherVolume(), G4ThreeVector(0., 0., 0.2*mm));

    holderSS1_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    holderSS2_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    holderSS3_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    holderSS4_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));

    beamline_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    beamlinevac_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    flange_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));

    cooling_log->SetVisAttributes(G4VisAttributes(G4Colour::White()));

    return nullptr;
}

void TargetChamberMg25::SetNewValue(G4UIcommand* command, G4String value)
{
    GeometryObject::SetNewValue(command, value);
}
