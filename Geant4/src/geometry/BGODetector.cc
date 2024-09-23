//BGODetector.cc

#include "geometry/BGODetector.hh"
#include "Analysis.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"

#include "G4SDManager.hh"
#include "G4Event.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDParticleFilter.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using CLHEP::um;
using CLHEP::mm;
using CLHEP::deg;

#include <memory>
using std::make_shared;


G4VPhysicalVolume* BGODetector::Construct()
{
    /// Dimensions
    //

    // inner aluminum tube
    const auto innerTubeRin    = 60.0 * mm * 0.5;
    const auto innerTubeRout   = innerTubeRin + 0.8 * mm;
    const auto innerTubeLength = (306 - 22) * mm;

    // steel mantle of the enclosure
    const auto mantleRin  = 240.0 * mm * 0.5;
    const auto mantleRout = 246.0 * mm * 0.5;
    const auto mantleLength = innerTubeLength;

    // outer cap
    const auto outerCapRin       =  60.0 * mm * 0.5;
    const auto outerCapRout      = 290.0 * mm * 0.5;
    const auto outerCapThickness =  11.0 * mm;

    // inner cap
    const auto innerCapRin       = 246.0 * mm * 0.5;
    const auto innerCapRout      = 290.0 * mm * 0.5;
    const auto innerCapThickness =  11.0 * mm;

    // crystal dimensions
    const auto shrinking = 0.2 * mm;
    const auto crystalLength     = 284.0 * mm;
    const auto crystalTrapShort  =  34.4 * mm - (shrinking * sin(60*deg));
    const auto crystalTrapLong   = 115.2 * mm - (shrinking * sin(60 * deg));
    const auto crystalTrapHeight =  70.0 * mm - (shrinking * cos(60 * deg));

    // aluminum BGO shielding
    const auto shieldingLength     = 284.0 * mm;
    const auto shieldingTrapShort  =  34.4 * mm;
    const auto shieldingTrapLong   = 115.2 * mm;
    const auto shieldingTrapHeight =  70.0 * mm;


    // cutout for PMTs in end caps of steel enclosure
    const auto pmtHoleRadius   = 0.5*53*mm;
    const auto pmtCircleRadius = 74*mm;

    const auto cuttingTolerance = 1*mm;

    /// Materials
    //
    auto* man = G4NistManager::Instance();

    auto* matBGO = man->FindOrBuildMaterial("G4_BGO");
    auto* matSteel  = G4Material::GetMaterial("Steel-316");
    auto* matAl = man->FindOrBuildMaterial("G4_Al");

    /// Solids
    //

    // Inner aluminum tube (surrounding borehole)
    auto* innerTube_tub =
        new G4Tubs(GetName() + "_innerTube",
                   innerTubeRin, innerTubeRout,
                   0.5*innerTubeLength,
                   0.0, 360.0*deg);

    // (Single) BGO crystal (trapezoid)
    auto crystal_trd =
        new G4Trd(GetName() + "_crystal",
                  0.5*crystalLength, 0.5*crystalLength,
                  0.5*crystalTrapLong, 0.5*crystalTrapShort, 0.5*crystalTrapHeight);

    // (Single) BGO shielding (trapezoid)
    auto shielding_trd =
        new G4Trd(GetName() + "_shielding",
                  0.5*shieldingLength, 0.5*shieldingLength,
                  0.5*shieldingTrapLong, 0.5*shieldingTrapShort, 0.5*shieldingTrapHeight);

    // Large outer steel cylinder
    auto* mantle_tub =
        new G4Tubs(GetName() + "_mantle",
                   mantleRin, mantleRout,
                   0.5*mantleLength,
                   0.0, 360.0*deg);

    // Outer part of end caps (down to borehole)
    auto* outerCap_tub =
        new G4Tubs(GetName() + "_outerCap",
                   outerCapRin, outerCapRout,
                   0.5*outerCapThickness,
                   0.0, 360.0*deg);

    // Inner part of end caps (down to cylinder)
    auto* innerCap_tub =
        new G4Tubs(GetName() + "_innerCap",
                   innerCapRin, innerCapRout,
                   0.5*innerCapThickness,
                   0.0, 360.0*deg);

    // Cut holes for PMTs into outer part of the endcaps
    // measurements (radius and position) are taken approximately from drawing, should be checked!
    G4VSolid *outerCap_solid = outerCap_tub;
    {
        auto* pmtHoleCutter =
            new G4Tubs(GetName() + "_pmtHoleCutter",
                       0, pmtHoleRadius,
                       0.5*outerCapThickness + cuttingTolerance,
                       0.0, 360.0*deg);

        for (int i = 0; i <= 5; i++)
        {
            G4RotationMatrix rotMat;
            rotMat.rotateZ(i*60*deg);
            G4ThreeVector posVec(pmtCircleRadius, 0.0*mm, 0.0*mm);
            posVec.rotateZ(i*60*deg);

            outerCap_solid = new G4SubtractionSolid(GetName() + "_endSteel1", outerCap_solid, pmtHoleCutter, &rotMat, posVec);
        }
    }

    /// Logical volumes
    //

    auto* innerTube_log =
        new G4LogicalVolume(innerTube_tub, matAl,
                            innerTube_tub->GetName());

    auto* mantle_log =
        new G4LogicalVolume(mantle_tub, matSteel,
                            mantle_tub->GetName());

    auto* outerCap_log =
        new G4LogicalVolume(outerCap_solid, matSteel,
                            outerCap_solid->GetName());

    auto* innerCap_log =
        new G4LogicalVolume(innerCap_tub, matSteel,
                            innerCap_tub->GetName());

    auto* crystal_log =
        new G4LogicalVolume(crystal_trd, matBGO,
                            crystal_trd->GetName());

    auto* shielding_log =
        new G4LogicalVolume(shielding_trd, matAl,
                            shielding_trd->GetName());

    /// Visualization attributes
    //
    innerTube_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    mantle_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    outerCap_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
    innerCap_log->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));

    crystal_log->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));
    shielding_log->SetVisAttributes(G4VisAttributes(G4Colour::Red()));


    /// Place as Physical Volumes
    //

    // Place six copies of BGO crystal, rotated by 60 degree each
    {
        G4RotationMatrix* bgoRot[6];
        G4ThreeVector     bgoVec[6];

        // Calculate positions and rotations
        for (int i = 0; i <= 5; i++)
        {
            bgoRot[i] = new G4RotationMatrix();
            bgoRot[i]->rotateY(90.0*deg);
            bgoRot[i]->rotateZ(i * 60*deg);
            bgoVec[i] = G4ThreeVector(-(30.8+70.0/2)*mm, 0.0*mm, 0.0*mm);
            bgoVec[i].rotateZ(i * 60*deg);
        }

        // Place volumes
	PlaceVolume(crystal_log, shielding_log,
		G4ThreeVector(0,0,0), G4RotationMatrix());
        for(G4int i = 0; i <= 5; i++)
	{
            PlaceVolume(shielding_log, GetMotherVolume(),
                        bgoVec[i], *bgoRot[i],
                        i); // copyNr
        }
    }


    // Enclosure
    {
        PlaceVolume(innerTube_log, GetMotherVolume(),
                    G4ThreeVector(0, 0, 0));

        PlaceVolume(mantle_log, GetMotherVolume(),
                    G4ThreeVector(0, 0, 0));

        // Upstream
        PlaceVolume(outerCap_log, GetMotherVolume(),
                    G4ThreeVector(0, 0, -0.5 * mantleLength - 0.5*outerCapThickness), G4RotationMatrix(),
                    0);

        // Downstream
        PlaceVolume(outerCap_log, GetMotherVolume(),
                    G4ThreeVector(0, 0,  0.5 * mantleLength + 0.5*outerCapThickness), G4RotationMatrix(),
                    1);

        // Upstream
        PlaceVolume(innerCap_log, GetMotherVolume(),
                    G4ThreeVector(0, 0, -0.5 * mantleLength + 0.5*innerCapThickness), G4RotationMatrix(),
                    0);

        // Downstream
        PlaceVolume(innerCap_log, GetMotherVolume(),
                    G4ThreeVector(0, 0,  0.5 * mantleLength - 0.5*innerCapThickness), G4RotationMatrix(),
                    1);
    }

    return nullptr;
}

void BGODetector::ConstructSDandField()
{
    // Energy deposition in the BGO crystals is tracked with a primitive scorer

    // attempts to delete this pointer in the destructor result in error,
    // apparently detectors are cleaned up by Geant4
    auto *det = new G4MultiFunctionalDetector("BGOCrystal");
    G4SDManager::GetSDMpointer()->AddNewDetector(det);

    auto *psEdep = new G4PSEnergyDeposit("Edep",1);
    det->RegisterPrimitive(psEdep);

    SetSensitiveDetector(GetName() + "_crystal", det);
}

void BGODetector::SetupOutput()
{
    auto *am = G4AnalysisManager::Instance();

    fTupleID = am->CreateNtuple("EdepBGO", "Energy Deposition in BGO detector");

    G4cout << "BGODetector::SetupOutput - " << fTupleID << G4endl;

    am->CreateNtupleDColumn("BGO1");
    am->CreateNtupleDColumn("BGO2");
    am->CreateNtupleDColumn("BGO3");
    am->CreateNtupleDColumn("BGO4");
    am->CreateNtupleDColumn("BGO5");
    am->CreateNtupleDColumn("BGO6");
    am->CreateNtupleDColumn("BGOsum");
    am->CreateNtupleDColumn("X");
    am->CreateNtupleDColumn("Y");
    am->CreateNtupleDColumn("Z");


    am->FinishNtuple();
}

void BGODetector::FillOutput(const G4Event *event)
{
    // Get hit collection ID
    if (fHCID == -1 )
    {
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID("BGOCrystal/Edep");
    }

    // Get sum values from hits collections
    //
    G4double Edep[6] = {0, 0, 0, 0, 0, 0};
    G4double EdepSum = 0.0;

    auto hitsCollection
        = static_cast<G4THitsMap<G4double>*>(
              event->GetHCofThisEvent()->GetHC(fHCID));

    for (auto it : *hitsCollection->GetMap())
    {
        Edep[it.first] = *(it.second);
        EdepSum += *(it.second);
    }

    // Fill ntuple
    auto *am = G4AnalysisManager::Instance();

    for (int i = 0; i < 6; i++)
    {
        am->FillNtupleDColumn(fTupleID, i, Edep[i] / CLHEP::MeV);
    }
    am->FillNtupleDColumn(fTupleID, 6, EdepSum / CLHEP::MeV);

    am->AddNtupleRow(fTupleID);
}

void BGODetector::FillPosition(G4ThreeVector &vec)
{

  auto *am = G4AnalysisManager::Instance();
  
  am->FillNtupleDColumn(fTupleID, 7, vec.getX());
  am->FillNtupleDColumn(fTupleID, 8, vec.getY());
  am->FillNtupleDColumn(fTupleID, 9, vec.getZ());

}
