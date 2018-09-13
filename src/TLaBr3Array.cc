#include "TLaBr3Array.hh"
#include "LaBrSD.hh"
#include "DetectorConstruction.hh"

#include "G4SDManager.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"

//#define DEBUG

#ifdef DEBUG
G4bool CHK_OVERLAP_LABR = true;
#else
G4bool CHK_OVERLAP_LABR = false;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TLaBr3Array::TLaBr3Array(const DetectorConstruction *aworld)
  : TDetector(aworld), fNumCrystal(8), fArraySD(0),
    fCrystalDiameter(8.89), fCrystalLength(20.32),
    fHousingThicknessFront(0.5), fHousingThicknessSide(2),
    fInsulationThicknessFront(0.05), fInsulationThicknessSide(0.05),
    fHousing(true), fTransparentHousing(false), fInsulation(true),
    fTransparentInsulation(false)
{
  fArrayPosition = G4ThreeVector(0,0,0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Spec()
void TLaBr3Array::Spec() const {
  G4cout << " ******** LaBr3Array:Spec() ********* " << G4endl
         << "Numer of crystals: " << fNumCrystal << G4endl
         << "Array Position: " << G4BestUnit(fArrayPosition,"Length") << G4endl
         << " ************************************ " << G4endl
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Build()
void TLaBr3Array::Build(){
  G4Material *LaBr3 = materialMgr->GetMaterial("LaBr3");
  G4Material *Al    = materialMgr->GetMaterial("Al");
  G4Material *MgO   = materialMgr->GetMaterial("MgO");

  G4ThreeVector Pos[8];
  G4double Psi[8], Theta[8], Phi[8];
  G4double crystalDiameter, crystalLength;

  G4int id;
  G4RotationMatrix Rot3D, Rot3D2;

  fGeomFile = fopen("./geometry/LaBr3/LaBr3GeometryIN.txt","r");
  if(fGeomFile==NULL){
    G4cout << "Cannot find the geometry file for LaBr3" << G4endl;
    return;
  }
  
  for(G4int i=0;!feof(fGeomFile) && i<fNumCrystal;i++){
    G4double x,y,z,psi,theta,phi; // position buffer
    char insulation[100];
    fscanf(fGeomFile,"%i %lf %lf %lf %lf %lf %lf %lf %lf %s",
	   &id,&x,&y,&z,&psi,&theta,&phi,&crystalDiameter,&crystalLength,
	   insulation);
    G4cout<<"Reading LaBr3 detectors "<<i<<": "<<x<<" "<<y<<" "<<z<<" "
          <<psi<<" "<<theta<<" "<<phi << " " << crystalDiameter << " "
	  << crystalLength << " " << insulation <<G4endl;
    Pos[i].set(x*cm, y*cm, z*cm);
    Psi[i] = psi;
    Theta[i] = theta;
    Phi[i] = phi;

    if(id==-1) break;
    if(id<1) continue;
  }


  // Crystal
  G4Tubs *solidLaBr3Crystal = new G4Tubs("LaBr3Crystal", // Name
					 0.0, // RMin
					 crystalDiameter*cm/2.0, // RMax
					 crystalLength*cm/2.0, // Dz
					 0.0*deg, // Angle min
					 360.0*deg); // Angle max

  G4LogicalVolume *logicLaBr3Crystal =
    new G4LogicalVolume(solidLaBr3Crystal, // solid
			LaBr3, // Material
			"logicLaBr3Crystal", // Name
			0,0,0);

  G4VisAttributes *visAttLaBr3Crystal =
    new G4VisAttributes(G4Colour(0.5,0.5,1.0));

  logicLaBr3Crystal->SetVisAttributes(visAttLaBr3Crystal);

  // Housing
  if(fHousing){
    G4double excessMaterialAtFront = 0.5; // cm
    G4double flangeLength=2.0, pmtLength=19.0;// cm

    G4Tubs *solidLaBr3HouseOUT = new G4Tubs("solidLaBr3HouseOUT",
					    0.0, //Rmin
					    (crystalDiameter+
					     2*fInsulationThicknessSide+
					     2*fHousingThicknessSide)*cm/2.,//RMax
					    (crystalLength+
					     2*fInsulationThicknessFront+
					     2*fHousingThicknessFront+
					     2.*excessMaterialAtFront)*cm/2.0,//Dz
					    0*deg, // Angle MIN
					    360*deg); // Angle MAX


    G4Tubs *solidLaBr3HouseIN = new G4Tubs("solidLaBr3HouseIn",
					   0.0,
					   (crystalDiameter+
					    2*fInsulationThicknessSide)*cm/2.0,
					   (crystalLength+
					    2*fInsulationThicknessFront+
					    2*fHousingThicknessFront+
					    2.0*excessMaterialAtFront+
					    pmtLength)*cm/2.0,
					   0*deg,
					   360*deg);

#ifdef DEBUG
    G4cout << "LaBr crystal Diameter is " << crystalDiameter
	 << ", lenght is " << crystalLength << G4endl;
    G4cout << "LaBr outer housing radius is: "
	   << (crystalDiameter+2*fInsulationThicknessSide+
	       2*fHousingThicknessSide)*cm/2.0
	   << " mm (should be 60.0 mm for Milano LaBr)" <<G4endl;
#endif

    //Back side: simplified structure
    G4Tubs *solidLaBr3HouseBackFlange = new G4Tubs("solidLaBr3HouseBackFlange",
						   0.0,
						   19.028*cm/2.0,
						   flangeLength*cm/2.0, 
						   0*deg,360*deg
						   );
    G4Tubs *solidLaBr3PMT = new G4Tubs("solidLaBr3PMT",
				       0.0,
				       12.0*cm/2.0,
				       pmtLength*cm/2.0,
				       0*deg,360*deg
				       );
    
    // Front cap
    G4Tubs *solidLaBr3FrontCap = new G4Tubs("LaBr3FrontCap",
					    0.0,
					    (crystalDiameter+
					     2*fInsulationThicknessSide)*cm/2.0,
					    fHousingThicknessFront*cm/2.0, 
					    0*deg,360*deg
					    );

    // Merge them
    G4UnionSolid *solidLaBr3HouseBackFlange2 =
      new G4UnionSolid("solidLaBr3HouseBackFlange2",
		       solidLaBr3HouseBackFlange,
		       solidLaBr3PMT,0,
		       G4ThreeVector(0.*cm,0.*cm,(flangeLength+pmtLength)*cm/2.)
		       );

    G4UnionSolid *solidLaBr3HouseBACKOUT =
      new G4UnionSolid("solidLaBr3HouseBACKOUT",
		       solidLaBr3HouseOUT,
		       solidLaBr3HouseBackFlange2,
		       0,
		       G4ThreeVector(0.,0.,(crystalLength+flangeLength)*cm/2.0)
		       );

    G4SubtractionSolid *solidLaBr3House1 =
      new G4SubtractionSolid("solidLaBr3HouseHouse1",
			     solidLaBr3HouseBACKOUT,
			     solidLaBr3HouseIN
			     );

    // Add Front cap
    G4UnionSolid *solidLaBr3House =
      new G4UnionSolid("solidLaBr3House",
		       solidLaBr3House1, solidLaBr3FrontCap,
		       0,
		       G4ThreeVector(0.,0.,
				     -(crystalLength+
				       2.0*fInsulationThicknessFront+
				       fHousingThicknessFront)*cm/2.0)
		       );


    // Logical Volume of the housing
    logicLaBr3House = new G4LogicalVolume(solidLaBr3House,
					  Al,
					  "logicLaBr3House",
					  0,0,0
					  );

    // vis attribute
    // red
    G4VisAttributes* visAttHouse =new G4VisAttributes(G4Colour(1.0,0.3,0.3)); 

    //set housing transparent
    if(fTransparentHousing) {visAttHouse->SetForceWireframe(true);} 
    logicLaBr3House->SetVisAttributes(visAttHouse);

  }// end of housing


  // Insulation
  if(fInsulation){
    G4Tubs *solidLaBr3InsulationOUT =
      new G4Tubs("solidLaBr3InsulationOUT",
		 0.0,
		 (crystalDiameter+2*fInsulationThicknessSide-0.0001)*cm/2.0,
		 (crystalLength+2*fInsulationThicknessFront-0.0001)*cm/2.0,
		 0*deg,
		 360*deg
		 );

    G4Tubs *solidLaBr3InsulationIN =
      new G4Tubs("solidLaBr3InsulationIN",
		 0.0,
		 crystalDiameter*cm/2.0,
		 crystalLength*cm/2.0,
		 0*deg,
		 360*deg
		 );

    G4SubtractionSolid *solidLaBr3Insulation = 
      new G4SubtractionSolid("solidLaBr3Insulation",
			     solidLaBr3InsulationOUT,
			     solidLaBr3InsulationIN
			     );
  
    logicLaBr3Insulation =  new G4LogicalVolume(solidLaBr3Insulation,
						MgO,
						"logicalLaBr3Insulation",
						0,0,0);

    G4VisAttributes* visAttInsulation =
      new G4VisAttributes(G4Colour(1.0,1.0,.0)); //yellow

    //set insulation transparent
    if(fTransparentInsulation) {visAttInsulation->SetForceWireframe(true);} 
    logicLaBr3Insulation->SetVisAttributes(visAttInsulation);
  } //insulation

  std::ofstream ofs;
  ofs.open("./geometry/LaBr3/LaBr3GeometryOUT.txt",
	   std::ios::out | std::ios::trunc);

  char tmp[100];
  for(G4int idet=0;idet<fNumCrystal;idet++){
    G4ThreeVector pos = Pos[idet] + fArrayPosition;
    G4cout << idet << ", position: " << pos << G4endl;
    
    ofs << std::setw(4) << std::left << idet + 1
        << std::setw(12) << std::setprecision(6) << std::left
        << pos.x() << " "
        << pos.y() << " "
        << pos.z() << G4endl;
    
    Rot3D.set(0, 0, 0);
    Rot3D.rotateX(Psi[idet]*degree);
    Rot3D.rotateY(Theta[idet]*degree);
    Rot3D.rotateZ(Phi[idet]*degree);

    // crystal
    sprintf(tmp,"physLaBr3Crystal%i",idet);
    G4VPhysicalVolume *physLaBr3Crystal = 
      new G4PVPlacement(G4Transform3D(Rot3D,pos), // 3D
                        tmp,// name
                        logicLaBr3Crystal,//logic
                        world->GetExpHall(),
                        FALSE,
                        idet,
                        CHK_OVERLAP_LABR
                        );

    // housing
    if(fHousing){
      sprintf(tmp,"physLaBr3House%i",idet);
      G4VPhysicalVolume *physLaBr3House = 
	new G4PVPlacement(G4Transform3D(Rot3D,pos),
			  tmp,
			  logicLaBr3House,
			  world->GetExpHall(),
			  FALSE,
			  idet+3000,
			  CHK_OVERLAP_LABR
			  );
    }
    
    // Insulation
    if(fInsulation){
      sprintf(tmp,"physLaBr3Insulation%i",idet);
      G4VPhysicalVolume *physLaBrInsulation = 
	new G4PVPlacement(G4Transform3D(Rot3D,pos),
			  tmp,
			  logicLaBr3Insulation,
			  world->GetExpHall(),
			  FALSE,
			  idet+4000,
			  CHK_OVERLAP_LABR
			  );
    }
  }
    
  ofs.close();

  G4SDManager *SDman = G4SDManager::GetSDMpointer();

  static G4int init3 = 1;
  if(init3){
    fArraySD = new LaBrSD("LaBrSD",
                          "LaBrCollection",
                          fNumCrystal);
    SDman->AddNewDetector(fArraySD);
    logicLaBr3Crystal->SetSensitiveDetector(fArraySD);
    init3 = 0;
  }

  static G4int init4 = 1;
  if(init4){
    fArrayRegion = new G4Region("LaBr");
    logicLaBr3Crystal->SetRegion(fArrayRegion);
    fArrayRegion->AddRootLogicalVolume(logicLaBr3Crystal);
    init4 = 0;
  }

}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TLaBr3Array::SetHousing(G4bool doHousing, G4double front, G4double side,
			     G4bool transparent){
  fHousingThicknessFront = front;
  fHousingThicknessSide = side;
  fTransparentHousing = transparent;

  G4cout << "The thickness of the LaBr3Array housing is: "
	 << fHousingThicknessFront << " "
	 <<fHousingThicknessSide << G4endl;
  
  if(fHousingThicknessFront<0 || fHousingThicknessSide<0){
    G4cout<<"HOUSING WILL BE OMITTED! " << G4endl;
    fHousing=false;
  }
  else{fHousing=doHousing;}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TLaBr3Array::SetInsulation(G4bool doInsulation,
				G4double front, G4double side,
				G4bool transparent){
  fInsulationThicknessFront = front;
  fInsulationThicknessSide = side;
  fTransparentInsulation = transparent;

  G4cout << "The thickness of the LaBr3Array insulation is: "
	 << fInsulationThicknessFront << " "
	 << fInsulationThicknessSide << G4endl;
  if(fInsulationThicknessFront<0 || fInsulationThicknessSide<0){
    G4cout << "INSULATION WILL BE OMITTED! "<< G4endl;
    fInsulation = false;
  }
  else{fInsulation = doInsulation;}
} 



/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float LaBr3Array::GetCrystalMeasuredEnergy(int a)  {
  if(fCrystalEnergy[a]==0) return 0.;
  if(fTypeOfEnergyResolution==0){return fCrystalEnergy[a];} //perfect resolution
  float dummy;
  if(fTypeOfEnergyResolution==1) dummy = (fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a]);
  else if(fTypeOfEnergyResolution==2) dummy = fEnergyResolution[0]*TMath::Power(fCrystalEnergy[a],fEnergyResolution[1]);
  else {
    if(fTypeOfEnergyResolution==3){
      dummy = TMath::Sqrt( fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a] + fEnergyResolution[2]*TMath::Power(fCrystalEnergy[a],2) );
      dummy/=2.35; //get sigma from fwhm
    }else{
    cout<<"Wrong resolution option '" << fTypeOfEnergyResolution << "' for LaBr3Array. Aborting program."<<endl; 
    abort();
    }
  }
  //Observed energy in the detector
  return G4RandGauss::shoot(fCrystalEnergy[a],dummy);
}
*/


