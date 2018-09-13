// TTestArray.cc
// 2012.11.8 Yasuhiro Togano

#include "TCATANA.hh"
#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"

#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
//#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"

#include <fstream>
#include "G4SystemOfUnits.hh"


//#define DEBUG

#ifdef DEBUG
G4bool CHK_OVERLAP = true;
G4int NumCrystalPerLayer = 1;
#else
G4bool CHK_OVERLAP = false;
G4int NumCrystalPerLayer = 20;
#endif


TCATANA::TCATANA(const DetectorConstruction *aworld)
  : TDetector(aworld), fArraySD(0),
    fAngCoveragePerCrystal(10*deg), fInRadius(25*cm), 
    fCrystalThickness(20*cm), fNumCrystal(1), 
    fMaterialName("CsI"), fBeamAcceptance(10*deg)
{
  fPosArray = G4ThreeVector(0,0,(-4635.14+138)*mm);
  G4cout<< ";;;;;;;;;;;;;;;;------------------;;;;;;;;;;;;;;;;;;;-----------------"<< G4endl
	<< "ZPosition = "<< fPosArray.z()/mm <<"mm"<<G4endl;
}

// Spec ======================================================================
void TCATANA::Spec() const {
  G4cout<< " ******** CATANA::Spec() ******** " << G4endl
	<< " Material: " << fMaterialName << G4endl
	<< " Inner Radius: " << G4BestUnit(fInRadius,"Length") << G4endl
	<< " Array center position: "<<G4BestUnit(fPosArray,"Length") << G4endl
	<< " Array center position: "<<fPosArray.z()/mm << G4endl
	<< " ************************************ " << G4endl
	<< G4endl;
}

// Build =====================================================================
void TCATANA::Build(){
  G4Material *crystal = materialMgr->GetMaterial(fMaterialName);
  G4Material *Al = materialMgr->GetMaterial("Al");

  // Length input: half length is a input
  G4double AngTheta = 10.5*deg;// theta coverage of this crystal
  G4double AngPhi = 18*deg; // phi coverage of this crystal
  G4double inrad = fInRadius; // init.macのsetInRad?
  G4double dx1 = inrad*std::sin(AngTheta/2.); // half length
  G4double pDz = 95/2.*mm; // crystal length
  G4double pTheta = 0*deg;
  G4double pPhi = 0*deg;
  G4double pDy1 = dx1;
  G4double pDx1 = inrad * std::cos(AngTheta) * std::tan(AngPhi/2.);
  G4double pDx2 = pDx1;
  G4double pAlp1 = 0*deg;
  G4double pDy2 = pDz*2. * tan(AngTheta/2.) + dx1;
  G4double pDx3 = (inrad+pDz*2/cos(AngTheta/2))*cos(AngTheta)* tan(AngPhi/2.);
  G4double pDx4 = pDx3;
  G4double pAlp2 = 0*deg;

  G4cout << "crystal outbox: pdy1: " << pDy1 << ", pdy2: "<< pDy2 
	 << ", pDx1: " << pDx1 << ", pDx3: "<< pDx3 << G4endl;

  G4Trap* fSolidHousingOUT1 = new G4Trap("SolidHousingOUT1",pDz, pTheta, pPhi,
					 pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					 pDx4, pAlp2);


  // 1mm thick Al housing
  pDz = 93./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 34.784/2.*mm;
  pDx1 = 60.606/2.*mm; 
  pDx2 = pDx1;
  pAlp1 = 0*deg;
  pDy2 = 51.875/2.*mm;
  pDx3 = 89.694/2.*mm;
  pDx4 = pDx3;
  pAlp2 = 0*deg;

  G4Trap* fSolidHousingIN1 = new G4Trap("SolidHousingIN1",pDz, pTheta, pPhi,
					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					pDx4, pAlp2);

  G4RotationMatrix rotHousing;
  rotHousing.set(0,0,0);
  G4SubtractionSolid* fSolidHousing1 = 
    new G4SubtractionSolid("SolidHousing1",fSolidHousingOUT1,// mother
			   fSolidHousingIN1, &rotHousing,
			   G4ThreeVector(0.,0.,0.));
  
  // Crystal
  pDz = 91./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 32.968/2.*mm;
  pDx1 = 58.919/2.*mm; 
  pDx2 = pDx1;
  pAlp1 = 0*deg;
  pDy2 = 49.692/2.*mm;
  pDx3 = 87.381/2.*mm;
  pDx4 = pDx3;
  pAlp2 = 0*deg;

  G4Trap* fSolidCrystal1 = new G4Trap("SolidCrystal1",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);

  G4LogicalVolume *fLogicHousing1 = new G4LogicalVolume(fSolidHousing1,
							Al,
							"logicHousing1",
							0,0,0);

  G4LogicalVolume *fLogicCrystal1 = new G4LogicalVolume(fSolidCrystal1,// solid
							crystal, // Material
							"logicCrystal1", // name
							0,0,0);

  G4int numdet = 0;  

  // geometry file output -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  std::ofstream ofs;
  ofs.open("./geometry/CATANA/CATANAGeometryOUT.txt",
	   std::ios::out | std::ios::trunc);
  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  G4ThreeVector crystalPos;
  
  G4double abs = inrad*cos(AngTheta/2) + pDz;
  for(G4int i=0;i<4;i++){
//    if(i==0){
//      G4double tabs = abs + 1.75*cm;
//      crystalPos.set(0.,tabs*cos(AngTheta/2.),-tabs*sin(AngTheta/2.));
//      crystalPos.rotateX(-AngTheta);
//    }
//    else if(i==1){
//      crystalPos.set(0.,abs*cos(AngTheta/2.)+2*mm,-abs*sin(AngTheta/2.)-0.2*mm);
//    }
    if(i==0 | i==1) continue;
    else if(i==2){
      crystalPos.set(0.,abs*cos(AngTheta/2.)+2*mm,abs*sin(AngTheta/2)+0.2*mm);
    }
    else if(i==3){
      G4double tabs = abs + 1.75*cm; // 15 cm version
      crystalPos.set(0.,tabs*cos(AngTheta/2),tabs*sin(AngTheta/2));
      crystalPos.rotateX(AngTheta);
    }
    for(G4int j=0;j<NumCrystalPerLayer;j++){
      G4RotationMatrix *rotCrystal = new G4RotationMatrix();
      rotCrystal->set(0.,0.,0.);
      //rotCrystal->rotateZ(-18*j*deg-9*deg);
      rotCrystal->rotateZ(-18*j*deg);
      rotCrystal->rotateX(105.75*deg- AngTheta*i);

#ifdef DEBUG
      G4cout << "ang: "<< 105.75-10.5*i << " deg, "<< numdet 
	     << "-th crystal: checking overlap" << G4endl;
      G4cout << " Volume of Crystal: "
	     << G4BestUnit(fSolidCrystal1->GetCubicVolume(),"Volume")
	     << G4endl;
#endif
      G4ThreeVector pos2 = crystalPos + fPosArray;
      char Pname[50];
      
      sprintf(Pname,"Housing%d",numdet);
      G4PVPlacement *fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
						     pos2,//pos
						     Pname,// Name
						     fLogicHousing1,
						     world->GetExpHall(),
						     FALSE, 
						     numdet+1000,
						     CHK_OVERLAP 
						     );

      
      sprintf(Pname,"Crystal%d",numdet);
      fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
				      pos2,//pos
				      Pname,// Name
				      fLogicCrystal1,
				      world->GetExpHall(),
				      FALSE, 
				      numdet,
				      CHK_OVERLAP 
				      );

      ofs << std::setw(4) <<  std::left << numdet + 1
	  << std::setw(12) << std::setprecision(6) << std::left
	  << pos2.x() << "  " 
	  << pos2.y() << "  " 
	  << pos2.z()
	  << G4endl;
    
      numdet++;
      crystalPos.rotateZ(18.*deg);
    }
  }

//=- crystal2 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Length input: half length is a input

  AngTheta = 10*deg;// theta coverage of this crystal
  AngPhi = 18*deg; // phi coverage of this crystal
  pDz = 105/2.*mm; // crystal length

  G4double detCY = 239.1244*mm;
  G4double detangle = 54.*deg;
  G4double faceY = detCY - pDz*sin(detangle);
  G4double rearY = detCY + pDz*sin(detangle);

  inrad = fInRadius;
  dx1 = inrad*std::sin(AngTheta/2.); // half length
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = dx1;
  pDy2 = pDz*2. * tan(AngTheta/2.) + dx1;
  pAlp1 = 0*deg;
  pAlp2 = 0*deg;
  pDx1 = (faceY - pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx2 = (faceY + pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx3 = (rearY - pDy2*cos(detangle))*tan(AngPhi/2.);
  pDx4 = (rearY + pDy2*cos(detangle))*tan(AngPhi/2.);

  G4cout << "crystal2 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
	 << G4endl
	 <<" pDy1: " << pDy1 << ", pDy2: "<< pDy2 << G4endl
	 <<" pDx1: "<< pDx1 << ", pDx2: " << pDx2 
	 << ", pDx3: " << pDx3 << ", pDx4: " << pDx4 <<G4endl;

  G4Trap* fSolidHousingOUT2 = new G4Trap("SolidHousingOUT2",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);



  // For 1mm thick Al housing
  pDz = 103./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 33.037/2.*mm;
  pDx1 = 57.217/2.*mm; 
  pDx2 = 63.368/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 51.050/2.*mm;
  pDx3 = 82.448/2.*mm;
  pDx4 = 91.955/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidHousingIN2 = new G4Trap("SolidHousingIN2",pDz, pTheta, pPhi,
					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					pDx4, pAlp2);

  rotHousing.set(0,0,0);
  G4SubtractionSolid* fSolidHousing2 = 
    new G4SubtractionSolid("SolidHousing2",fSolidHousingOUT2,// mother
			   fSolidHousingIN2, &rotHousing,
			   G4ThreeVector(0.,0.,0.));


  // Crystal
  pDz = 101./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 31.212/2.*mm;
  pDx1 = 55.387/2.*mm; 
  pDx2 = 61.199/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 48.903/2.*mm;
  pDx3 = 80.650/2.*mm;
  pDx4 = 89.756/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidCrystal2 = new G4Trap("SolidCrystal2",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);

  G4LogicalVolume *fLogicHousing2 = new G4LogicalVolume(fSolidHousing2,
							Al,
							"logicHousing2",
							0,0,0);

  G4LogicalVolume *fLogicCrystal2 = new G4LogicalVolume(fSolidCrystal2,// solid
							crystal, // Material
							"logicCrystal2", // name
							0,0,0);

  for(G4int i=0;i<2;i++){
    if(i==0){
      //crystalPos.set(0.,313.44*mm + 1*mm,140.39*mm + 1*mm); // 15 cm ver
      //crystalPos.set(0., 241.4184*mm+0.1*mm, 108.2998*mm);
      crystalPos.set(0., 251.479*mm, 120.9738*mm);
    }
    else if(i==1){
      //crystalPos.set(0.,298.77*mm+25*mm, 192.89*mm+13*mm); // 15 cm ver
      crystalPos.set(0., 239.1244*mm, 166.1174*mm);
    }
    for(G4int j=0;j<NumCrystalPerLayer;j++){
      G4RotationMatrix *rotCrystal = new G4RotationMatrix();
      rotCrystal->set(0.,0.,0.);
      rotCrystal->rotateZ(-18*j*deg);
      rotCrystal->rotateX(64*deg- AngTheta*i);

#ifdef DEBUG
      G4cout << "ang: "<< 64-10*i << " deg, "<< numdet 
	     << "-th crystal: checking overlap" << G4endl;
      G4cout << " Volume of Crystal: "
	     << G4BestUnit(fSolidCrystal2->GetCubicVolume(),"Volume")
	     << G4endl;
#endif
      G4ThreeVector pos2 = crystalPos + fPosArray;

      char Pname[50];
      
      sprintf(Pname,"Housing%d",numdet);
      G4PVPlacement *fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
						     pos2,
						     Pname,// Name
						     fLogicHousing2,
						     world->GetExpHall(),
						     FALSE, 
						     numdet+1000,
						     CHK_OVERLAP 
						     );

      
      sprintf(Pname,"Crystal%d",numdet);
      fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
				      pos2,
				      Pname,// Name
				      fLogicCrystal2,
				      world->GetExpHall(),
				      FALSE, 
				      numdet,
				      CHK_OVERLAP 
				      );


      ofs << std::setw(4) <<  std::left << numdet + 1
	  << std::setw(12) << std::setprecision(6) << std::left
	  << pos2.x() << "  " 
	  << pos2.y() << "  " 
	  << pos2.z()
	  << G4endl;
    
      numdet++;
      crystalPos.rotateZ(18.*deg);
    }
  }

//=- crystal3 =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  // Length input: half length is a input

  AngTheta = 11*deg;// theta coverage of this crystal
  AngPhi = 18*deg; // phi coverage of this crystal
  pDz = 125/2.*mm; // crystal length

  detCY = 215.5226*mm;
  detangle = 43.5*deg;
  faceY = detCY - pDz*sin(detangle);
  rearY = detCY + pDz*sin(detangle);

  inrad = fInRadius;
  dx1 = inrad*std::sin(AngTheta/2.); // half length
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = dx1;
  pDy2 = pDz*2. * tan(AngTheta/2.) + dx1;
  pAlp1 = 0*deg;
  pAlp2 = 0*deg;
  pDx1 = (faceY - pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx2 = (faceY + pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx3 = (rearY - pDy2*cos(detangle))*tan(AngPhi/2.);
  pDx4 = (rearY + pDy2*cos(detangle))*tan(AngPhi/2.);

  G4cout << "crystal3 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
	 << G4endl
	 <<" pDy1: " << pDy1 << ", pDy2: "<< pDy2 << G4endl
	 <<" pDx1: "<< pDx1 << ", pDx2: " << pDx2 
	 << ", pDx3: " << pDx3 << ", pDx4: " << pDx4 <<G4endl;

  G4Trap* fSolidHousingOUT3 = new G4Trap("SolidHousingOUT3",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);


  // For 1mm thick Al housing
  pDz = 123./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 36.531/2.*mm;
  pDx1 = 48.446/2.*mm; 
  pDx2 = 56.840/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 60.218/2.*mm;
  pDx3 = 72.981/2.*mm;
  pDx4 = 86.817/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidHousingIN3 = new G4Trap("SolidHousingIN3",pDz, pTheta, pPhi,
					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					pDx4, pAlp2);

  rotHousing.set(0,0,0);
  G4SubtractionSolid* fSolidHousing3 = 
    new G4SubtractionSolid("SolidHousing3",fSolidHousingOUT3,// mother
			   fSolidHousingIN3, &rotHousing,
			   G4ThreeVector(0.,0.,0.));


  // Crystal
  pDz = 121./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 34.724/2.*mm;
  pDx1 = 46.654/2.*mm; 
  pDx2 = 54.632/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 58.025/2.*mm;
  pDx3 = 71.233/2.*mm;
  pDx4 = 84.565/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidCrystal3 = new G4Trap("SolidCrystal3",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);

  G4LogicalVolume *fLogicHousing3 = new G4LogicalVolume(fSolidHousing3,
							Al,
							"logicHousing3",
							0,0,0);

  G4LogicalVolume *fLogicCrystal3 = new G4LogicalVolume(fSolidCrystal3,// solid
							crystal, // Material
							"logicCrystal3", // name
							0,0,0);

  crystalPos.set(0.,215.5226*mm, 208.8918*mm);
  for(G4int j=0;j<NumCrystalPerLayer;j++){
    G4RotationMatrix *rotCrystal = new G4RotationMatrix();
    rotCrystal->set(0.,0.,0.);
    rotCrystal->rotateZ(-18*j*deg);
    rotCrystal->rotateX(detangle);

#ifdef DEBUG
    G4cout << "ang: "<< detangle/degree << " deg, "<< numdet 
	   << "-th crystal: checking overlap" << G4endl;
    G4cout << " Volume of Crystal: "
	   << G4BestUnit(fSolidCrystal3->GetCubicVolume(),"Volume")
	   << G4endl;
#endif
    G4ThreeVector pos2 = crystalPos + fPosArray;

    char Pname[50];
    
    sprintf(Pname,"Housing%d",numdet);
    G4PVPlacement *fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
						   pos2,
						   Pname,// Name
						   fLogicHousing3,
						   world->GetExpHall(),
						   FALSE, 
						   numdet+1000,
						   CHK_OVERLAP 
						   );

    sprintf(Pname,"Crystal%d",numdet);
    fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
				    pos2,
				    Pname,// Name
				    fLogicCrystal3,
				    world->GetExpHall(),
				    FALSE, 
				    numdet,
				    CHK_OVERLAP 
				    );
    

    ofs << std::setw(4) <<  std::left << numdet + 1
	<< std::setw(12) << std::setprecision(6) << std::left
	<< pos2.x() << "  " 
	<< pos2.y() << "  " 
	<< pos2.z()
	<< G4endl;
    
    numdet++;
    crystalPos.rotateZ(18.*deg);
  }
    //  }


  // crystal4 =================================================================
  AngTheta = 11*deg;// theta coverage of this crystal
  AngPhi = 18*deg; // phi coverage of this crystal
  pDz = 140/2.*mm; // crystal length

  //detCY = 171.9715*mm;
  detCY = 180.0394*mm;
  detangle = 32.5*deg;
  faceY = detCY - pDz*sin(detangle);
  rearY = detCY + pDz*sin(detangle);

  inrad = fInRadius;
  dx1 = inrad*std::sin(AngTheta/2.); // half length
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = dx1;
  pDy2 = pDz*2. * tan(AngTheta/2.) + dx1;
  pAlp1 = 0*deg;
  pAlp2 = 0*deg;
  pDx1 = (faceY - pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx2 = (faceY + pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx3 = (rearY - pDy2*cos(detangle))*tan(AngPhi/2.);
  pDx4 = (rearY + pDy2*cos(detangle))*tan(AngPhi/2.);

  G4cout << "crystal4 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
	 << G4endl
	 <<" pDy1: " << pDy1 << ", pDy2: "<< pDy2 << G4endl
	 <<" pDx1: "<< pDx1 << ", pDx2: " << pDx2 
	 << ", pDx3: " << pDx3 << ", pDx4: " << pDx4 <<G4endl;

  G4Trap* fSolidHousingOUT4 = new G4Trap("SolidHousingOUT4",pDz, pTheta, pPhi,
					 pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					 pDx4, pAlp2);

  // For 1mm thick Al housing
  pDz = 138./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 36.531/2.*mm;
  pDx1 = 38.237/2.*mm; 
  pDx2 = 47.997/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 63.082/2.*mm;
  pDx3 = 58.515/2.*mm;
  pDx4 = 75.368/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidHousingIN4 = new G4Trap("SolidHousingIN4",pDz, pTheta, pPhi,
					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					pDx4, pAlp2);

  rotHousing.set(0,0,0);
  G4SubtractionSolid* fSolidHousing4 = 
    new G4SubtractionSolid("SolidHousing4",fSolidHousingOUT4,// mother
			   fSolidHousingIN4, &rotHousing,
			   G4ThreeVector(0.,0.,0.));


  // Crystal
  pDz = 136./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 34.627/2.*mm;
  pDx1 = 36.478/2.*mm; 
  pDx2 = 45.730/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 61.010/2.*mm;
  pDx3 = 56.782/2.*mm;
  pDx4 = 73.082/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidCrystal4 = new G4Trap("SolidCrystal4",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);

  G4LogicalVolume *fLogicHousing4 = new G4LogicalVolume(fSolidHousing4,
							Al,
							"logicHousing4",
							0,0,0);

  G4LogicalVolume *fLogicCrystal4 = new G4LogicalVolume(fSolidCrystal4,// solid
							crystal, // Material
							"logicCrystal4", // name
							0,0,0);

  //crystalPos.set(0., 171.9715*mm, 266.4462*mm);
  crystalPos.set(0., 180.0394*mm, 247.0881*mm);
  //for(G4int j=0;j<0;j++){
  for(G4int j=0;j<NumCrystalPerLayer;j++){
    G4RotationMatrix *rotCrystal = new G4RotationMatrix();
    rotCrystal->set(0.,0.,0.);
    rotCrystal->rotateZ(-18*j*deg);
    rotCrystal->rotateX(detangle);

#ifdef DEBUG
    G4cout << "ang: "<< detangle/degree << " deg, "<< numdet 
	   << "-th crystal: checking overlap" << G4endl;
    G4cout << " Volume of Crystal: "
	   << G4BestUnit(fSolidCrystal4->GetCubicVolume(),"Volume")
	   << G4endl;
#endif
    G4ThreeVector pos2 = crystalPos + fPosArray;

    char Pname[50];
    sprintf(Pname,"Housing%d",numdet);
    G4PVPlacement *fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
						   pos2,
						   Pname,// Name
						   fLogicHousing4,
						   world->GetExpHall(),
						   FALSE, 
						   numdet+1000,
						   CHK_OVERLAP 
						   );
    sprintf(Pname,"Crystal%d",numdet);
    fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
				    pos2,
				    Pname,// Name
				    fLogicCrystal4,
				    world->GetExpHall(),
				    FALSE, 
				    numdet,
				    CHK_OVERLAP 
				    );

    ofs << std::setw(4) <<  std::left << numdet + 1
	<< std::setw(12) << std::setprecision(6) << std::left
	<< pos2.x() << "  " 
	<< pos2.y() << "  " 
	<< pos2.z()
	<< G4endl;
    
    numdet++;
    crystalPos.rotateZ(18.*deg);
  }

  // crystal5 =================================================================
  AngTheta = 11*deg;// theta coverage of this crystal
  AngPhi = 18*deg; // phi coverage of this crystal
  pDz = 150/2.*mm; // crystal length

  detCY = 134.9135*mm;
  detangle = 21.50*deg;
  faceY = detCY - pDz*sin(detangle);
  rearY = detCY + pDz*sin(detangle);

  inrad = fInRadius;
  dx1 = inrad*std::sin(AngTheta/2.); // half length
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = dx1;
  pDy2 = pDz*2. * tan(AngTheta/2.) + dx1;
  pAlp1 = 0*deg;
  pAlp2 = 0*deg;
  pDx1 = (faceY - pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx2 = (faceY + pDy1*cos(detangle))*tan(AngPhi/2.);
  pDx3 = (rearY - pDy2*cos(detangle))*tan(AngPhi/2.);
  pDx4 = (rearY + pDy2*cos(detangle))*tan(AngPhi/2.);

  G4cout << "crystal5 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
	 << G4endl
	 <<" pDy1: " << pDy1 << ", pDy2: "<< pDy2 << G4endl
	 <<" pDx1: "<< pDx1 << ", pDx2: " << pDx2 
	 << ", pDx3: " << pDx3 << ", pDx4: " << pDx4 <<G4endl;

  G4Trap* fSolidHousingOUT5 = new G4Trap("SolidHousingOUT5",pDz, pTheta, pPhi,
					 pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					 pDx4, pAlp2);

  // For 1mm thick Al housing
  pDz = 148./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 36.531/2.*mm;
  pDx1 = 26.646/2.*mm; 
  pDx2 = 37.412/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 65.032/2.*mm;
  pDx3 = 39.860/2.*mm;
  pDx4 = 59.027/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidHousingIN5 = new G4Trap("SolidHousingIN5",pDz, pTheta, pPhi,
					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
					pDx4, pAlp2);

  rotHousing.set(0,0,0);
  G4SubtractionSolid* fSolidHousing5 = 
    new G4SubtractionSolid("SolidHousing5",fSolidHousingOUT5,// mother
			   fSolidHousingIN5, &rotHousing,
			   G4ThreeVector(0.,0.,0.));


  // Crystal
  pDz = 146./2.*mm;
  pTheta = 0*deg;
  pPhi = 0*deg;
  pDy1 = 34.724/2.*mm;
  pDx1 = 24.912/2.*mm; 
  pDx2 = 35.146/2.*mm;
  pAlp1 = 0*deg;
  pDy2 = 62.840/2.*mm;
  pDx3 = 38.183/2.*mm;
  pDx4 = 56.704/2.*mm;
  pAlp2 = 0*deg;

  G4Trap* fSolidCrystal5 = new G4Trap("SolidCrystal5",pDz, pTheta, pPhi,
				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
				      pDx4, pAlp2);

  G4LogicalVolume *fLogicHousing5 = new G4LogicalVolume(fSolidHousing5,
							Al,
							"logicHousing5",
							0,0,0);

  G4LogicalVolume *fLogicCrystal5 = new G4LogicalVolume(fSolidCrystal5,// solid
							crystal, // Material
							"logicCrystal5", // name
							0,0,0);

  //crystalPos.set(0., 130.7253*mm, 291.2625*mm);
  crystalPos.set(0., 134.9135*mm, 275.1572*mm);
  for(G4int j=0;j<NumCrystalPerLayer;j++){
    //for(G4int j=0;j<0;j++){
    G4RotationMatrix *rotCrystal = new G4RotationMatrix();
    rotCrystal->set(0.,0.,0.);
    rotCrystal->rotateZ(-18*j*deg);
    rotCrystal->rotateX(detangle);

#ifdef DEBUG
    G4cout << "ang: "<< detangle/degree << " deg, "<< numdet 
	   << "-th crystal: checking overlap" << G4endl;
    G4cout << " Volume of Crystal: "
	   << G4BestUnit(fSolidCrystal5->GetCubicVolume(),"Volume")
	   << G4endl;
#endif
    G4ThreeVector pos2 = crystalPos + fPosArray;
    char Pname[50];
    sprintf(Pname,"Housing%d",numdet);
    G4PVPlacement *fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
						   pos2,
						   Pname,// Name
						   fLogicHousing5,
						   world->GetExpHall(),
						   FALSE, 
						   numdet+1000,
						   CHK_OVERLAP 
						   );
    sprintf(Pname,"Crystal%d",numdet);
    fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
				    pos2,
				    Pname,// Name
				    fLogicCrystal5,
				    world->GetExpHall(),
				    FALSE, 
				    numdet,
				    CHK_OVERLAP 
				    );

    ofs << std::setw(4) <<  std::left << numdet + 1
	<< std::setw(12) << std::setprecision(6) << std::left
	<< pos2.x() << "  " 
	<< pos2.y() << "  " 
	<< pos2.z()
	<< G4endl;
    
    numdet++;
    crystalPos.rotateZ(18.*deg);
  }

  // crystal6 =================================================================
  // AngTheta = 11*deg;// theta coverage of this crystal
  // AngPhi = 18*deg; // phi coverage of this crystal
  // pDz = 150/2.*mm; // crystal length

  // detCY = 83.41*mm;
  // detangle = 10.5*deg;
  // faceY = detCY - pDz*sin(detangle);
  // rearY = detCY + pDz*sin(detangle);

  // inrad = fInRadius;
  // dx1 = inrad*std::sin(AngTheta/2.); // half length
  // pTheta = 0*deg;
  // pPhi = 0*deg;
  // pDy1 = dx1;
  // pDy2 = pDz*2. * tan(AngTheta/2.) + dx1;
  // pAlp1 = 0*deg;
  // pAlp2 = 0*deg;
  // pDx1 = (faceY - pDy1*cos(detangle))*tan(AngPhi/2.);
  // pDx2 = (faceY + pDy1*cos(detangle))*tan(AngPhi/2.);
  // pDx3 = (rearY - pDy2*cos(detangle))*tan(AngPhi/2.);
  // pDx4 = (rearY + pDy2*cos(detangle))*tan(AngPhi/2.);

  // G4cout << "crystal6 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" 
  // 	 << G4endl
  // 	 <<" pDy1: " << pDy1 << ", pDy2: "<< pDy2 << G4endl
  // 	 <<" pDx1: "<< pDx1 << ", pDx2: " << pDx2 
  // 	 << ", pDx3: " << pDx3 << ", pDx4: " << pDx4 <<G4endl;

  // G4Trap* fSolidHousingOUT6 = new G4Trap("SolidHousingOUT6",pDz, pTheta, pPhi,
  // 					 pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
  // 					 pDx4, pAlp2);

  // // For 1mm thick Al housing
  // pDz = 148./2.*mm;
  // pTheta = 0*deg;
  // pPhi = 0*deg;
  // pDy1 = 36.531/2.*mm;
  // pDx1 = 14.403/2.*mm; 
  // pDx2 = 25.781/2.*mm;
  // pAlp1 = 0*deg;
  // pDy2 = 65.073/2.*mm;
  // pDx3 = 18.611/2.*mm;
  // pDx4 = 38.879/2.*mm;
  // pAlp2 = 0*deg;

  // G4Trap* fSolidHousingIN6 = new G4Trap("SolidHousingIN6",pDz, pTheta, pPhi,
  // 					pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
  // 					pDx4, pAlp2);

  // rotHousing.set(0,0,0);
  // G4SubtractionSolid* fSolidHousing6 = 
  //   new G4SubtractionSolid("SolidHousing6",fSolidHousingOUT6,// mother
  // 			   fSolidHousingIN6, &rotHousing,
  // 			   G4ThreeVector(0.,0.,0.));


  // // Crystal
  // pDz = 146./2.*mm;
  // pTheta = 0*deg;
  // pPhi = 0*deg;
  // pDy1 = 34.724/2.*mm;
  // pDx1 = 12.685/2.*mm; 
  // pDx2 = 23.500/2.*mm;
  // pAlp1 = 0*deg;
  // pDy2 = 62.840/2.*mm;
  // pDx3 = 16.965/2.*mm;
  // pDx4 = 36.537/2.*mm;
  // pAlp2 = 0*deg;

  // G4Trap* fSolidCrystal6 = new G4Trap("SolidCrystal6",pDz, pTheta, pPhi,
  // 				      pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3,
  // 				      pDx4, pAlp2);

  // G4LogicalVolume *fLogicHousing6 = new G4LogicalVolume(fSolidHousing6,
  // 							Al,
  // 							"logicHousing6",
  // 							0,0,0);

  // G4LogicalVolume *fLogicCrystal6 = new G4LogicalVolume(fSolidCrystal6,// solid
  // 							crystal, // Material
  // 							"logicCrystal6", // name
  // 							0,0,0);


//   crystalPos.set(0.,83.41*mm, 289.6388*mm);
//   //for(G4int j=0;j<0;j++){
//   for(G4int j=0;j<NumCrystalPerLayer;j++){
//     G4RotationMatrix *rotCrystal = new G4RotationMatrix();
//     rotCrystal->set(0.,0.,0.);
//     rotCrystal->rotateZ(-18*j*deg);
//     rotCrystal->rotateX(detangle);

// #ifdef DEBUG
//     G4cout << "ang: "<< detangle/degree << " deg, "<< numdet 
// 	   << "-th crystal: checking overlap" << G4endl;
//     G4cout << " Volume of Crystal: "
// 	   << G4BestUnit(fSolidCrystal6->GetCubicVolume(),"Volume")
// 	   << G4endl;
// #endif
//     G4ThreeVector pos2 = crystalPos + fPosArray;
//     char Pname[50];
//     sprintf(Pname,"Housing%d",numdet);
//     G4PVPlacement *fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
// 						   pos2,
// 						   Pname,// Name
// 						   fLogicHousing6,
// 						   world->GetExpHall(),
// 						   FALSE, 
// 						   numdet+1000,
// 						   CHK_OVERLAP 
// 						   );
//     sprintf(Pname,"Crystal%d",numdet);
//     fPhysiArray = new G4PVPlacement(rotCrystal, // rotation
// 				    pos2,
// 				    Pname,// Name
// 				    fLogicCrystal6,
// 				    world->GetExpHall(),
// 				    FALSE, 
// 				    numdet,
// 				    CHK_OVERLAP 
// 				    );

//     ofs << std::setw(4) <<  std::left << numdet + 1
// 	<< std::setw(12) << std::setprecision(6) << std::left
// 	<< pos2.x() << "  " 
// 	<< pos2.y() << "  " 
// 	<< pos2.z()
// 	<< G4endl;
    
//     numdet++;
//     crystalPos.rotateZ(18.*deg);
//  }

  
  G4cout << "Total Number of Crystal:" << numdet << G4endl;
  
  fNumCrystal = numdet;

  ofs.close();

  //Visualization
  //G4VisAttributes *fAttArray0 = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  ////fAttArray0->SetVisibility(false);
  //fLogicCrystal0->SetVisAttributes(fAttArray0);
  G4VisAttributes *fAttArray1 = new G4VisAttributes(G4Colour(0.1,0.5,1.0));
  //fAttArray1->SetVisibility(false);
  fLogicHousing1->SetVisAttributes(fAttArray1);
  G4VisAttributes *fAttArray2 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  //fAttArray2->SetVisibility(false);
  fLogicHousing2->SetVisAttributes(fAttArray2);
  G4VisAttributes *fAttArray3 = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  //fAttArray3->SetVisibility(false);
  fLogicHousing3->SetVisAttributes(fAttArray3);
  G4VisAttributes *fAttArray4 = new G4VisAttributes(G4Colour(1.0,0.5,0.0));
  //fAttArray4->SetVisibility(false);
  fLogicCrystal4->SetVisAttributes(fAttArray4);
  G4VisAttributes *fAttArray5 = new G4VisAttributes(G4Colour(1.0,0.75,0.0));
  //fAttArray5->SetVisibility(false);
  fLogicCrystal5->SetVisAttributes(fAttArray5);
  //  G4VisAttributes *fAttArray6 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  //fAttArray6->SetVisibility(false);
  //  fLogicCrystal6->SetVisAttributes(fAttArray6);

    //Region
  static int init=1;
  if(init){
    fArrayRegion = new G4Region("Array");
    //fLogicCrystal0->SetRegion(fArrayRegion);
    fLogicCrystal1->SetRegion(fArrayRegion);
    fLogicCrystal2->SetRegion(fArrayRegion);
    fLogicCrystal3->SetRegion(fArrayRegion);
    fLogicCrystal4->SetRegion(fArrayRegion);
    fLogicCrystal5->SetRegion(fArrayRegion);
    //    fLogicCrystal6->SetRegion(fArrayRegion);
    //fArrayRegion->AddRootLogicalVolume(fLogicCrystal0);
    fArrayRegion->AddRootLogicalVolume(fLogicCrystal1);
    fArrayRegion->AddRootLogicalVolume(fLogicCrystal2);
    fArrayRegion->AddRootLogicalVolume(fLogicCrystal3);
    fArrayRegion->AddRootLogicalVolume(fLogicCrystal4);
    fArrayRegion->AddRootLogicalVolume(fLogicCrystal5);
    //    fArrayRegion->AddRootLogicalVolume(fLogicCrystal6);
    init = 0;
  }

  //Sensitive Detector
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  static int init2=1;
  if(init2){
    fArraySD = new CalorimeterSD("arraySD",
				 "arrayCollection",
				 numdet);
    SDman->AddNewDetector(fArraySD);

    fLogicCrystal1->SetSensitiveDetector(fArraySD);
    fLogicCrystal2->SetSensitiveDetector(fArraySD);
    fLogicCrystal3->SetSensitiveDetector(fArraySD);
    fLogicCrystal4->SetSensitiveDetector(fArraySD);
    fLogicCrystal5->SetSensitiveDetector(fArraySD);
    //    fLogicCrystal6->SetSensitiveDetector(fArraySD);
    init2 = 0;
  }


}
