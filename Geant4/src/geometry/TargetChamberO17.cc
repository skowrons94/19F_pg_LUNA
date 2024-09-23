#include "geometry/TargetChamberO17.hh"

#include <list>

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


TargetChamberO17::TargetChamberO17() : GeometryObject("TargetChamberO17")
{
}

TargetChamberO17::~TargetChamberO17()
{
}

G4VPhysicalVolume *TargetChamberO17::Construct()
{

	/// Dimensions
	//

	//CF100 flange upstream

	const auto lengthCF100 = 18 * mm;
	const auto rInnerCF100 = 0.5 * 37 * mm; // Assumed everything is a cylinder with an inner radius of 37 mm
	const auto rOuterCF100 = 0.5 * 114 * mm;

	/// CF40 flange

	const auto lengthCF40 = 13 * mm;
	const auto rInnerCF40 = 0.5 * 37 * mm;
	const auto rOuterCF40 = 0.5 * 70 * mm;

	///  Tube

	const auto lengthTube = 234 * mm;
	const auto rInnerTube = 0.5 * 37 * mm;
	const auto rOuterTube = 0.5 * 40 * mm;

	/// Threaded counter flange for target holder initialization

	const auto lengthCounterFlange = 13 * mm; // was 12.7 * mm
	const auto rInnerCounterFlange = 0.5 * 38 * mm;
	const auto rOuterCounterFlange = 0.5 * 52 * mm;

	///  Cold finger
	//// Assuming 1 mm thickness

	const auto lengthColdFinger = lengthTube + lengthCounterFlange + lengthCF40 + lengthCF100;
	const auto rInnerColdFinger = 0.5 * 24 * mm;
	const auto rOuterColdFinger = 0.5 * 26 * mm;


	const auto phiStart = 0*deg;
	const auto phiTotal = 360*deg;

	/// Chamber itself

	G4int chamberNZPlanes = 4;
	const G4double chamberZPlanes[] =
	{
		lengthCF100+lengthCF40,
		lengthCF100+lengthCF40+lengthTube,
		lengthCF100+lengthCF40+lengthTube,
		lengthCF100+lengthCF40+lengthTube+lengthCounterFlange
	};

	const G4double chamberRInner[] =
	{
		rInnerTube,
		rInnerTube,
		rInnerCounterFlange,
		rInnerCounterFlange
	};

	const G4double chamberROuter[] =
	{
		rOuterTube,
		rOuterTube,
		rOuterCounterFlange,
		rOuterCounterFlange
	};

	/// Materials
	//

	G4NistManager* man = G4NistManager::Instance();

	G4Material* matSteel = G4Material::GetMaterial("Steel-316");
	G4Material* matAl = man->FindOrBuildMaterial("G4_Al");
	G4Material* matCopper = man->FindOrBuildMaterial("G4_Cu");

	auto matChamber = matAl;

	/// Solids
	//

	auto solidCF100
		= new G4Tubs("CF100",
			rInnerCF100, rOuterCF100,
			0.5*lengthCF100,
			phiStart,phiTotal);

	auto solidCF40
		= new G4Tubs("CF40",
			rInnerCF40, rOuterCF40,
			0.5*lengthCF40,
			phiStart,phiTotal);

	G4Tubs* solidColdFinger
		= new G4Tubs("ColdFinger",
				rInnerColdFinger,rOuterColdFinger,
				0.5 * lengthColdFinger,
				phiStart,phiTotal);

	G4Polycone* solidTargetChamber = new G4Polycone("TargetChamber",phiStart,phiTotal,
			chamberNZPlanes,chamberZPlanes,chamberRInner,chamberROuter);

	/// Logical volumes
	//

	G4LogicalVolume* logicCF100 = new G4LogicalVolume(solidCF100, matSteel,
			solidCF100->GetName());

	G4LogicalVolume* logicCF40 = new G4LogicalVolume(solidCF40, matSteel,
			solidCF40->GetName());

	G4LogicalVolume* logicColdFinger = new G4LogicalVolume(solidColdFinger, matCopper,
			solidColdFinger->GetName());

	G4LogicalVolume* logicTargetChamber = new G4LogicalVolume(solidTargetChamber, matChamber,
			solidTargetChamber->GetName());

	/// Visualization attributes
	//

	logicCF100->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
	logicCF40->SetVisAttributes(G4VisAttributes(G4Colour::Green()));
	logicColdFinger->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));
	logicTargetChamber->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));

	/// Place as Physical Volumes
	//
	const G4double targetLength = 0.25*mm;

	PlaceVolume(logicCF100, GetMotherVolume(),
			-G4ThreeVector(0.0,0.0, 0.5*lengthCF100 + lengthCF40 + lengthTube + lengthCounterFlange-targetLength));

	PlaceVolume(logicCF40, GetMotherVolume(),
			-G4ThreeVector(0.0,0.0, 0.5*lengthCF40 + lengthTube + lengthCounterFlange-targetLength));

	PlaceVolume(logicColdFinger, GetMotherVolume(),
			-G4ThreeVector(0.0,0.0, 0.5 * lengthColdFinger + 13 * mm ));

	PlaceVolume(logicTargetChamber, GetMotherVolume(),
			-G4ThreeVector(0.0,0.0, lengthCF100 + lengthCF40 + lengthTube + lengthCounterFlange-targetLength));

	return nullptr;
}

void TargetChamberO17::SetNewValue(G4UIcommand* command, G4String value)
{
    GeometryObject::SetNewValue(command, value);
}
