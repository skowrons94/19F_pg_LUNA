#include "geometry/BeamLine.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"

#include "G4Tubs.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using CLHEP::mm;
using CLHEP::um;
using CLHEP::deg;


G4VPhysicalVolume *BeamLine::Construct()
{
    /// Dimensions
    //
    const auto beamLineLength = 125.0*mm;
    const auto flangeWidth    =  12.5*mm; // ??

    const auto tolerance = 10*um;


    /// Materials
    //
    auto* man = G4NistManager::Instance();

    auto* matVacuum = man->FindOrBuildMaterial("G4_Galactic");
    auto* matSteel  = G4Material::GetMaterial("Steel-316");


    /// Solids
    //
    auto* beamline_tub =
        new G4Tubs(GetName() + "_Tube",
                   23.5*mm, 25.0*mm,               // Rin, Rout
                   0.5*beamLineLength - tolerance, // half-length
                   0.0, 360.*deg);                 // phi0, d(phi)

    auto* beamlinevac_tub =
        new G4Tubs(GetName() + "_beamlinevac2",
                   0, 23.5*mm - tolerance,         // Rin, Rout
                   0.5*beamLineLength,             // half-length
                   0.0, 360.*deg);                 // phi0, d(phi)

    auto* flange_tub =
        new G4Tubs(GetName() + "_flange2",
                   25*mm + tolerance, 50.0*mm,     // Rin, Rout
                   0.5*flangeWidth,                // half-length
                   0.0, 360.*deg);                 // phi0, d(phi)


    /// Logical Volumes
    //
    auto* beamline_log =
        new G4LogicalVolume(beamline_tub, matSteel,
                            beamline_tub->GetName());

    auto* beamlinevac_log =
        new G4LogicalVolume(beamlinevac_tub, matVacuum,
                            beamlinevac_tub->GetName());

    auto* flange_log =
        new G4LogicalVolume(flange_tub, matSteel,
                            flange_tub->GetName());


    /// Visualization Attributes
    //
    beamline_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    beamlinevac_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    flange_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));


    /// Place as Physical Volumes
    //
    PlaceVolume(beamline_log,      GetMotherVolume(), -G4ThreeVector(0.,0.,0.5*beamLineLength + 280.*mm));
    PlaceVolume(beamlinevac_log,   GetMotherVolume(), -G4ThreeVector(0.,0.,0.5*beamLineLength + 280.*mm));

    PlaceVolume(flange_log,   GetMotherVolume(),   -G4ThreeVector(0.,0.,0.5*flangeWidth + 280*mm));

    return nullptr;
}
