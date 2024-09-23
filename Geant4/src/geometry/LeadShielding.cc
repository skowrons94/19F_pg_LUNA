#include "geometry/LeadShielding.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4SubtractionSolid.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using CLHEP::mm;
using CLHEP::cm;
using CLHEP::deg;

#include <iostream>
using std::cout;
using std::endl;

#include <memory>
using std::make_shared;

LeadShielding::LeadShielding() :
    GeometryObject("LeadShielding")
{
    fSetupCmd = make_shared<G4UIcmdWithAString>("/Geometry/LeadShielding/setup", this);
    fSetupCmd->SetGuidance("Set setup of the lead shielding.");
    fSetupCmd->SetParameterName("mode", false);
    fSetupCmd->SetCandidates("BGO Germanium0deg Germanium55deg");
    fSetupCmd->SetToBeBroadcasted(false);

    f55DegreePositionCmd = make_shared<G4UIcmdWithADoubleAndUnit>("/Geometry/LeadShielding/55degreePosition", this);
    f55DegreePositionCmd->SetGuidance("Set setup of the lead shielding.");
    f55DegreePositionCmd->SetParameterName("position", false);
    f55DegreePositionCmd->SetToBeBroadcasted(false);
}

void LeadShielding::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if(command == fSetupCmd.get())
    {
        fSetupChosen = true;

        if (newValue == "BGO")
        {
            fSetup = setupBGO;
        }
        else if (newValue == "Germanium0deg")
        {
            fSetup = setupGermanium0deg;
        }
        else if (newValue == "Germanium55deg")
        {
            fSetup = setupGermanium55deg;
        }
        else
        {
            G4cout << "Unknown LeadShielding setup '" << newValue << "'." << G4endl;
            fSetupChosen = false;
        }
    }
    else if (command == f55DegreePositionCmd.get())
    {
        f55degreePosition = f55DegreePositionCmd->GetNewDoubleValue(newValue);
    }
    else
    {
        GeometryObject::SetNewValue(command, newValue);
    }
}


G4VPhysicalVolume *LeadShielding::Construct()
{

    // Reminder:
    //  +x is towards to gas target beamline
    //  +y is upwards
    //  +z is beam direction

    if (!fSetupChosen)
    {
        cout << "No setup chosen for LeadShielding. Please use /Geometry/LeadShielding/setup command." << endl;
        exit(99);
    }

    G4NistManager* man = G4NistManager::Instance();
    G4Material* matLead = man->FindOrBuildMaterial("G4_Pb");

    G4RotationMatrix noRotation;

    const double extraSize = 1*cm; // extra dimension to be used for boolean solids to avoid sharing surfaces

    // start with a solid block of lead that contains upstream and downstream part

    const double upStreamLength = 305.5*mm;
    const double downStreamLength = 468*mm;

    const double leadBlockLengthX = 50*cm;
    const double leadBlockLengthY = (22.5+25)*cm;
    const double leadBlockLengthZ = upStreamLength + downStreamLength;

    const double separationZ = 0.5*leadBlockLengthZ - downStreamLength; // z-position of "gap" between the two parts, measured from the upstream end of the lead

    G4VSolid *leadShieldingSolid =
        new G4Box(GetName() + "_body",
                  0.5*leadBlockLengthX,
                  0.5*leadBlockLengthY,
                  0.5*leadBlockLengthZ);

    // first drill out the holes along the beam-axis
    const double innerChamberRadius = 147.0*mm;

    double beamCanalZ[6];
    beamCanalZ[0] = 0; // upstream, outside
    beamCanalZ[1] = beamCanalZ[0] + extraSize + (305.5-155.5)*mm; // downstream end of canal for flange
    beamCanalZ[2] = beamCanalZ[1]; // as above
    beamCanalZ[3] = beamCanalZ[2] + (155.5 + 318)*mm; // downstream end of large inner chamber
    beamCanalZ[4] = beamCanalZ[3]; // as above
    beamCanalZ[5] = beamCanalZ[4] + 150*mm + extraSize; // downstream, outside

    // relative to the center of the big lead block
    const double beamCanalY = -12.5*mm;

    double beamCanalInnerRadii[6] = {0, 0, 0, 0, 0, 0};
    double beamCanalOuterRadii[6] = {57.5*mm, 57.5*mm, innerChamberRadius*mm, innerChamberRadius*mm, 25*mm, 25*mm};
    G4Polycone *beamCanalSolid = new G4Polycone(GetName() + "_beamCanal", 0, 360.0*deg, 6, beamCanalZ, beamCanalInnerRadii, beamCanalOuterRadii);


    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, beamCanalSolid,
            G4Transform3D(noRotation,
                          G4ThreeVector(0, beamCanalY, -0.5*leadBlockLengthZ-extraSize)));

    // vector to the center of the inner chamber
    const G4ThreeVector innerChamberCenter(0, beamCanalY, -0.5*leadBlockLengthZ+0.5*(beamCanalZ[2]+beamCanalZ[3])-extraSize);

    // we need to place many flat cuts now. for this all we need is a sufficiently large box to use and re-use.
    const double cutBoxLengthX = 220*mm; // large enough for the largest required cut
    const double cutBoxLengthY = 220*mm;
    const double cutBoxLengthZ = leadBlockLengthZ + 2*extraSize;
    G4Box *cutBoxSolid = new G4Box(GetName() + "_cutBox", 0.5*cutBoxLengthX, 0.5*cutBoxLengthY, 0.5*cutBoxLengthZ);

    G4RotationMatrix rotate45z;
    rotate45z.rotateZ(45*deg);

    // "top" = +Y
    // "right" = +Z (downstream)
    // "back" = +X (towards other beamline)

    // top back
    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, cutBoxSolid,
            G4Transform3D(rotate45z,
                          G4ThreeVector(             105*mm + 0.5*145*mm + 0.5*cutBoxLengthX/sqrt(2),
                                  beamCanalY + 105*mm + 0.5*145*mm + 0.5*cutBoxLengthY/sqrt(2),
                                  0)));
    // parallel bottom back
    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, cutBoxSolid,
            G4Transform3D(noRotation,
                          G4ThreeVector(117*mm + 0.5*cutBoxLengthX,
                                        beamCanalY - 192*mm - 0.5*cutBoxLengthY,
                                        0)));

    // parallel bottom front
    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, cutBoxSolid,
            G4Transform3D(noRotation,
                          G4ThreeVector(-117*mm - 0.5*cutBoxLengthX,
                                        beamCanalY - 192*mm - 0.5*cutBoxLengthY,
                                        0)));


    // sloped bottom back
    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, cutBoxSolid,
            G4Transform3D(rotate45z,
                          G4ThreeVector(             165*mm + 0.5*85*mm + 0.5*cutBoxLengthX/sqrt(2),
                                  beamCanalY - 107*mm - 0.5*85*mm - 0.5*cutBoxLengthY/sqrt(2),
                                  0)));

    // now the parts that are specific to the up-/downstream sections

    // upstream top front
    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, cutBoxSolid,
            G4Transform3D(rotate45z,
                          G4ThreeVector(           - 150*mm - 0.5*100*mm - 0.5*cutBoxLengthX/sqrt(2),
                                  beamCanalY + 150*mm + 0.5*100*mm + 0.5*cutBoxLengthY/sqrt(2),
                                  separationZ - 0.5*cutBoxLengthZ)));

    // downstream top front
    leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, cutBoxSolid,
            G4Transform3D(rotate45z,
                          G4ThreeVector(           - 210*mm - 0.5*40*mm - 0.5*cutBoxLengthX/sqrt(2),
                                  beamCanalY + 210*mm + 0.5*40*mm + 0.5*cutBoxLengthY/sqrt(2),
                                  0))); // technically separationZ + 0.5*cutBoxLengthZ))), but upstream has larger cut-out anyway

    // now we drill the small canals for the water cooling
    const double verticalCoolingTubeZ = -0.5*leadBlockLengthZ+upStreamLength+5*mm;

    G4Box *drillWaterCoolingSquareSolid =
        new G4Box(GetName() + "_drillWaterCoolingSquare", 5*mm, 0.5*leadBlockLengthY+extraSize, 2.5*mm);
    leadShieldingSolid =
        new G4SubtractionSolid(GetName() + "_body",
                               leadShieldingSolid,
                               drillWaterCoolingSquareSolid,
                               &noRotation, G4ThreeVector(0, 0, verticalCoolingTubeZ-2.5*mm));
    G4RotationMatrix rotate90x = noRotation;
    rotate90x.rotateX(90*deg);
    G4Tubs *drillWaterCoolingTubeSolid =
        new G4Tubs(GetName() + "_drillWaterCoolingTube", 0*mm, 5*mm, 0.5*leadBlockLengthY+extraSize, 0.0, 360.0*deg);
    leadShieldingSolid =
        new G4SubtractionSolid(GetName() + "_body",
                               leadShieldingSolid,
                               drillWaterCoolingTubeSolid,
                               G4Transform3D(rotate90x, G4ThreeVector(0, 0, verticalCoolingTubeZ)));


    // for the Germanium55deg setup one will need to drill the exit here.
    const double drill55degLength = 150*mm; // suitable for drilling both holes
    G4Tubs *drill55deg = new G4Tubs(GetName() + "_drill55deg", 0, 65*mm, drill55degLength, 0, 360.0*deg);
    G4RotationMatrix rotate305y = noRotation;
    rotate305y.rotateY(305*deg);

    if (fSetup == setupGermanium55deg && false)
    {
        // the 55 degree hole in the main lead piece is open
        leadShieldingSolid = new G4SubtractionSolid(GetName() + "_body", leadShieldingSolid, drill55deg,
                G4Transform3D(rotate305y,
                              G4ThreeVector((-250+0.5*105)*mm,
                                            beamCanalY,
                                            separationZ + 183.6*mm - 0.5*105*mm/tan(55*deg))));
    }

    // move detector so that the target can be at 0, 0, 0
    G4ThreeVector shift = -innerChamberCenter;

    if (fSetup == setupGermanium55deg && false)
    {
        shift += G4ThreeVector(10.1575*mm, 0, 76.25*mm - 11.8*mm*cos(55*deg));
    }
    else if (fSetup == setupBGO)
    {
        shift += G4ThreeVector(0, 0, 76.25*mm);
    }


    // place the shielding
    G4LogicalVolume* blockLog =
        new G4LogicalVolume(leadShieldingSolid, matLead,
                            leadShieldingSolid->GetName());
    blockLog->SetVisAttributes(G4VisAttributes(G4Colour(0.8,0.8,0.8)));

    PlaceVolume(blockLog, GetMotherVolume(), shift);


    // for the Germanium setup we have an additional inset, with an extra exit at 55°
    if (fSetup == setupGermanium0deg || fSetup == setupGermanium55deg)
    {
        G4VSolid *insetSolid = new G4Tubs(GetName() + "_inset1", 100.0*mm, 145*mm, 0.5*471.5*mm, 0, 360.0*deg);
        // the inset is slightly smaller than the surrounding shielding, hence no overlap (extraSize) is required, but it's a separate volume
        insetSolid = new G4SubtractionSolid(GetName() + "_inset2", insetSolid, drill55deg,
                                            G4Transform3D(rotate305y,
                                                    G4ThreeVector(-145*mm,
                                                            0,
                                                            -0.5*471.5*mm+264.5*mm)));

        // drilling the grooves for the water cooling
        G4Tubs *insetCutoutTubeSolid = new G4Tubs(GetName() + "_insetCutoutTube", 0*mm, 5*mm, 150*mm, 0, 360.0*deg);
        insetSolid = new G4SubtractionSolid(GetName() + "_inset3", insetSolid, insetCutoutTubeSolid, G4Transform3D(rotate90x, G4ThreeVector(0,0,-76.25*mm)));

        G4Box *insetCutoutBoxSolid = new G4Box(GetName() + "_inserCutoutBox", 5*mm, 150*mm, 159.5*mm+extraSize);
        insetSolid = new G4SubtractionSolid(GetName() + "_inset4", insetSolid, insetCutoutBoxSolid, &noRotation, G4ThreeVector(0,0, -0.5*471.5*mm-insetCutoutBoxSolid->GetZHalfLength() + 159.5*mm));

        G4LogicalVolume* insetLog =
            new G4LogicalVolume(insetSolid, matLead,
                                insetSolid->GetName());
        insetLog->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));

        PlaceVolume(insetLog, GetMotherVolume(), shift + innerChamberCenter);
    }

    if (f55degreePosition >= 0 && false)
    {
        G4VSolid *sidePieceSolid =
            new G4CutTubs(GetName() + "_sidePieceStep1",
                          0.5*130.0*mm,
                          0.5*420.0*mm,
                          0.5*167.0*mm,
                          0.0, 360.0*deg,
                          G4ThreeVector(0,0,-1),
                          G4ThreeVector(-sin(35*deg),0,cos(35*deg)));

        G4LogicalVolume *sidePieceLog =
            new G4LogicalVolume(sidePieceSolid, matLead, sidePieceSolid->GetName());
        sidePieceLog->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

        G4RotationMatrix rotate125y = noRotation;
        rotate125y.rotateY(125*deg);
        PlaceVolume(sidePieceLog,
                    GetMotherVolume(),
                    shift + f55degreePosition * G4ThreeVector(-sin(55*deg),0,cos(55*deg)) + G4ThreeVector(-25*cm-16.7*cm*0.5*cos(35*deg), beamCanalY, -0.5*leadBlockLengthZ+305.5*mm+183.6*mm+0.5*167.0*mm*sin(35*deg)), rotate125y);
    }

    return nullptr;
}


