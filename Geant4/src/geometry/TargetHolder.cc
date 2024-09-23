#include "geometry/TargetHolder.hh"

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

/* #include <stdexcept> */
/* using std::runtime_error; */

TargetHolder::TargetHolder() : GeometryObject("TargetHolder")
{
}

TargetHolder::~TargetHolder()
{
}

G4VPhysicalVolume *TargetHolder::Construct()
{
	/// Materials
	//
	G4NistManager* man = G4NistManager::Instance();

	G4Material* matAl = man->FindOrBuildMaterial("G4_Al");
	G4Material* matBrass = man->FindOrBuildMaterial("G4_BRASS");
	G4Material* matWater = man->FindOrBuildMaterial("G4_WATER");
	G4Material* matTa = man->FindOrBuildMaterial("G4_Ta");
	G4Material* matC = man->FindOrBuildMaterial("G4_GRAPHITE");
	G4Material* matSteel = G4Material::GetMaterial("Steel-316");
	G4Material* matAlloy = man->FindOrBuildMaterial("Alloy");

	auto matHolder = matAl;

	/// Dimensions
	//

	G4double lengthHolderLower;
	G4double rInnerHolderLower;
	G4double rOuterHolderLower;

	G4double lengthHolderMidLow;
	G4double rInnerHolderMidLow;
	G4double rOuterHolderMidLow;

	G4double lengthHolderMid;
	G4double rInnerHolderMid;
	G4double rOuterHolderMid;

	G4double lengthHolderUpper;
	G4double rInnerHolderUpper;
	G4double rOuterHolderUpper;

	G4double lengthHolderTop;
	G4double rInnerHolderTop;
	G4double rOuterHolderTop;

	G4double lengthTarget;
	G4double rInnerTarget;
	G4double rOuterTarget;

	G4double lengthWater;
	G4double rInnerWater;
	G4double rOuterWater;

	G4double lengthWaterMid;
	G4double rInnerWaterMid;
	G4double rOuterWaterMid;

	G4double lengthWaterUpper;
	G4double rInnerWaterUpper;
	G4double rOuterWaterUpper;

	G4double phiStart;
	G4double phiTotal;

	// Graphite Target
	const G4double lengthTargetGraph = -6.0 * mm;
	const G4double rInnerTargetGraph = 0.5 * 0 * mm;
	const G4double rOuterTargetGraph = 0.5 * 25 * mm;

	if(matHolder == matAl)
	{
	        lengthHolderLower = 7 * mm; // was 10 * mm
		rInnerHolderLower = 0.5 * 52 * mm + 1 * um;
		rOuterHolderLower = 0.5 * 58 * mm; // was 59 * mm

		lengthHolderMidLow = 3 * mm;      // new
		rInnerHolderMidLow = 0.5 * 54 * mm; // new
		rOuterHolderMidLow = 0.5 * 58 * mm; // new

		lengthHolderMid = 1.7 * mm;
		rInnerHolderMid = 0.5 * 23.2 * mm;
		rOuterHolderMid = 0.5 * 58 * mm; // was 59 * mm

		lengthHolderUpper = 7.5 * mm; // was 9.2 * mm
		rInnerHolderUpper = 0.5 * 23.2 * mm;
		rOuterHolderUpper = 0.5 * 40.5 * mm; // was 38 * mm

		lengthHolderTop = 3.5 * mm; // was (11-9.2) * mm
		rInnerHolderTop = 0.5 * 0 * mm;
		rOuterHolderTop = 0.5 * 40.5 * mm; // was 38 * mm

		lengthTarget = 0.25 * mm;
		rInnerTarget = 0.5 * 0 * mm;
		rOuterTarget = 0.5 * 41 * mm; // was 40.5 * mm

		lengthWater = 9.2 * mm - 100 * um; // was (9.2+1.7) * mm
		rInnerWater = 0.5 * 0 * mm;
		rOuterWater = 0.5 * 23.2 * mm - 100 * um;

		phiStart = 0*deg;
		phiTotal = 360*deg;
	}
	if(matHolder == matBrass)
	{
		lengthHolderLower = 10 * mm;
		rInnerHolderLower = 0.5 * 50.3 * mm;
		rOuterHolderLower = 0.5 * 57 * mm;

		lengthHolderMid = (19.7-14.0) * mm;
		rInnerHolderMid = 0.5 * 29.5 * mm;
		rOuterHolderMid = 0.5 * 57 * mm;

		lengthHolderUpper = (14.0-5.3) * mm;
		rInnerHolderUpper = 0.5 * 23.2 * mm;
		rOuterHolderUpper = 0.5 * 41.2 * mm;

		lengthHolderTop = 5.3 * mm;
		rInnerHolderTop = 0.5 * 0 * mm;
		rOuterHolderTop = 0.5 * 41.2 * mm;

		 lengthTarget = 0.25 * mm;
		 rInnerTarget = 0.5 * 0 * mm;
		 rOuterTarget = 0.5 * 30 * mm;


		 lengthWaterMid = lengthHolderMid -100 * um;
		 rInnerWaterMid = 0.5 * 0 * mm;
		 rOuterWaterMid = 0.5 * 29.5 * mm -100 * um;

		 lengthWaterUpper = lengthHolderUpper -100 * um;
		 rInnerWaterUpper = 0.5 * 0 * mm;
		 rOuterWaterUpper = 0.5 * 23.2 * mm -100 * um;

		 phiStart = 0*deg;
		 phiTotal = 360*deg;
	}

	/// ORing
	const G4double lengthORing = 1 * mm;
	const G4double rInnerORing = 0.5 * 27 * mm;
	const G4double rOuterORing = 0.5 * 40 * mm;

	/// Target Holder
	const  G4int  holderNZPlanes = 10;
	const  G4double holderZPlanes[] =
	{
		0,
		lengthHolderLower,
		lengthHolderLower,
		lengthHolderLower+lengthHolderMidLow,
		lengthHolderLower+lengthHolderMidLow,
		lengthHolderLower+lengthHolderMid+lengthHolderMidLow,
		lengthHolderLower+lengthHolderMid+lengthHolderMidLow,
		lengthHolderLower+lengthHolderMid+lengthHolderUpper+lengthHolderMidLow,
		lengthHolderLower+lengthHolderMid+lengthHolderUpper+lengthHolderMidLow,
		lengthHolderLower+lengthHolderMid+lengthHolderUpper+lengthHolderTop+lengthHolderMidLow
	};

	const G4double holderRInner[] =
	{
		rInnerHolderLower,
		rInnerHolderLower,
		rInnerHolderMidLow,
		rInnerHolderMidLow,
		rInnerHolderMid,
		rInnerHolderMid,
		rInnerHolderUpper,
		rInnerHolderUpper,
		rInnerHolderTop,
		rInnerHolderTop
	};

	const G4double holderROuter[] =
	{
		rOuterHolderLower,
		rOuterHolderLower,
		rOuterHolderMidLow,
		rOuterHolderMidLow,
		rOuterHolderMid,
		rOuterHolderMid,
		rOuterHolderUpper,
		rOuterHolderUpper,
		rOuterHolderTop,
		rOuterHolderTop
	};

	/// Graphite Target
	const G4int targetGraphiteNZPlanes = 2;
	const G4double targetGraphiteZPlanes[] =
	{
		0,
		lengthTargetGraph
	};
	const G4double targetGraphiteRInner[] = {rInnerTargetGraph,rInnerTargetGraph};
	const G4double targetGraphiteROuter[] = {rOuterTargetGraph,rOuterTargetGraph};

	/// Target
	const G4int targetNZPlanes = 2;
	const G4double targetZPlanes[] =
	{
		0,
		lengthTarget
	};
	const G4double targetRInner[] = {rInnerTarget,rInnerTarget};
	const G4double targetROuter[] = {rOuterTarget,rOuterTarget};

	/// Volumes
	//

	/// Solids
	G4Polycone* solidHolder = new G4Polycone("Holder",phiStart,phiTotal,holderNZPlanes,
			holderZPlanes,holderRInner,holderROuter);

	G4Polycone* solidTarget = new G4Polycone("Target",phiStart,phiTotal,targetNZPlanes,
			targetZPlanes,targetRInner,targetROuter);

	G4Polycone* solidTargetGraphite = new G4Polycone("TargetGraphite",
							 phiStart,phiTotal,
							 targetGraphiteNZPlanes,
							 targetGraphiteZPlanes,
							 targetGraphiteRInner,
							 targetGraphiteROuter);

	G4Tubs* solidORing = new G4Tubs("ORing", rInnerORing,rOuterORing,
				0.5 * lengthORing, phiStart, phiTotal);

	/// Logical
	G4LogicalVolume* logicHolder = new G4LogicalVolume(solidHolder, matAlloy,
			solidHolder->GetName());

	G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, matTa,
			solidTarget->GetName());

	G4LogicalVolume* logicTargetGraphite = new G4LogicalVolume(solidTargetGraphite,
								   matC,
								   solidTargetGraphite->GetName());

	G4LogicalVolume* logicORing = new G4LogicalVolume(solidORing, matSteel,
			solidORing->GetName());

	/// Water
	if(matHolder == matAl){
		const G4int waterNZPlanes = 2;
		const G4double waterZPlanes[] =
			{
				lengthHolderMidLow,
				lengthHolderMidLow + lengthWater
			};
		const G4double waterRInner[] = {rInnerWater, rInnerWater};
		const G4double waterROuter[] = {rOuterWater, rOuterWater};

		G4Polycone* solidWater = new G4Polycone("Water",phiStart,phiTotal,waterNZPlanes,
			waterZPlanes,waterRInner,waterROuter);

		G4LogicalVolume* logicWater = new G4LogicalVolume(solidWater, matWater,
			solidWater->GetName());

		logicWater->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

		PlaceVolume(logicWater, GetMotherVolume(),
			G4ThreeVector(0.,0.,lengthTarget));

	} else if(matHolder == matBrass){
		const G4int waterNZPlanes = 4;
		const G4double waterZPlanes[] =
			{
				0,
				lengthWaterMid,
				lengthWaterMid,
				lengthWaterMid + lengthWaterUpper
			};
		const G4double waterRInner[] = {rInnerWaterMid,rInnerWaterMid,
			rInnerWaterUpper, rInnerWaterUpper};
		const G4double waterROuter[] = {rOuterWaterMid, rOuterWaterMid,
			rOuterWaterUpper, rOuterWaterUpper};

		G4Polycone* solidWater = new G4Polycone("Water",phiStart,phiTotal,waterNZPlanes,
			waterZPlanes,waterRInner,waterROuter);

		G4LogicalVolume* logicWater = new G4LogicalVolume(solidWater, matWater,
			solidWater->GetName());

		logicWater->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

		PlaceVolume(logicWater, GetMotherVolume(),
			G4ThreeVector(0.,0.,lengthTarget));
	}

	/// Visualization attributes
	logicHolder->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
	logicTarget->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
	logicORing->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));

	/// Place as Physical Volumes
	PlaceVolume(logicHolder, GetMotherVolume(),
			-G4ThreeVector(0.,0., holderZPlanes[1]-lengthTarget));

	PlaceVolume(logicTarget, GetMotherVolume(),
			-G4ThreeVector(0.,0.,-lengthHolderMidLow));

	//	PlaceVolume(logicTargetGraphite, GetMotherVolume(),
	//	      -G4ThreeVector(0.,0.,0.));

	PlaceVolume(logicORing, GetMotherVolume(),
			-G4ThreeVector(0.,0., -lengthHolderMidLow + 0.5 * lengthORing + 1 * um));

	return nullptr;
}

void TargetHolder::SetNewValue(G4UIcommand* command, G4String value)
{
    GeometryObject::SetNewValue(command, value);
}
