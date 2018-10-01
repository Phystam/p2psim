// DetectorConstruction for Gamma detector 
// 2012.10.15 Yasuhiro Togano

#include "DetectorConstruction.hh"
#include "MaterialManager.hh"

#include "G4SDManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "simconst.hh"
#include "G4Trd.hh"
#include "G4SystemOfUnits.hh"
#include "TString.h"
#include "SiSD.hh"
// define the size of the experimental hall
static const G4ThreeVector SIZE_EXPHALL(5.*m, 5.*m, 10*m);

// class descriptions
DetectorConstruction::DetectorConstruction()
  : solidExpHall(0), logicExpHall(0), physiExpHall(0),
    size_ExpHall(SIZE_EXPHALL), fCATANA(0),
    fLaBr3Array(0), fInRadius(25*cm), fDetMaterialName("CsI")
{
  fArrayPos.set(0,0,0);
  fLaBrPos.set(0,0,0);
  fDetMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete fDetMessenger;
  delete fCATANA;
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  //Cleanup old geometory for the DetectorMessenger ----------------------
  //Without cleanup, the Detector Control based on UI does not work
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  MaterialManager *materialMgr = MaterialManager::GetPointer();

  // Experimental Hall
  solidExpHall = new G4Box("ExpHall_solid", size_ExpHall.x()/2.,
			   size_ExpHall.y()/2., size_ExpHall.z()/2.);

  G4Material *Air = materialMgr -> GetMaterial("Vacuum");

  logicExpHall = new G4LogicalVolume(solidExpHall, Air, "ExpHall_logic");

  // Visualization
  G4VisAttributes *expHallVisAtt = new G4VisAttributes(FALSE, G4Color(1,1,1));
  logicExpHall -> SetVisAttributes(expHallVisAtt);

  // placement of Experimental Hall
  physiExpHall = new G4PVPlacement(0,               // no lotation
				   G4ThreeVector(), // at (0,0,0)
				   logicExpHall,    // logicalVolume
				   "ExpHall_phys",  // name
				   0,               // mother volume
				   false,           // no boolean operation
				   0);              // copy number
  
  // Detectors ======================================================

    
  fCATANA = new TCATANA(this);
  fCATANA->SetInnerRadius(fInRadius);
  fCATANA->SetPositionVector(fArrayPos);
  fCATANA->SetMaterialName(fDetMaterialName);
  fCATANA->SetPositionVector(fArrayPos);
  fCATANA->Spec();
  fCATANA->Build();
  
  //LaBr3s ===========================================================

  //  fLaBr3Array = new TLaBr3Array(this);
  //  fLaBr3Array->SetArrayPosition(fLaBrPos);
  //  fLaBr3Array->Spec();
  //  fLaBr3Array->Build();

  // Si trackers ====================================================
  G4double tgt_center_z= TARGETPOS;

  G4Material* Si = materialMgr->GetMaterial("Si");
  G4double SiX = 50.*mm;
  G4double SiY = 50.*mm;
  G4double SiZ = 50.*um;
  const G4int numz = 2;
  G4double zplane[]={0.*mm, 40.*mm};
  G4double rinner[]={100.*mm,40.*mm};
  G4double router[]={100.3*mm,40.3*mm};
  //  G4Box* SiBox = new G4Box("SiBox", SiX, SiY, SiZ);
  G4Polyhedra* SiPoly = new G4Polyhedra("SiPoly",0,60.*deg,1,numz,zplane,rinner,router);

  // G4Trd* SiPoly = new G4Trd("SiPoly",-10*mm,10*mm,-20*mm,20*mm,300*um);

  G4LogicalVolume* SiLV = new G4LogicalVolume(SiPoly, Si, "SiLV", 0, 0, 0);

  G4int numdet =0;
  for(int isi=0;isi<6;isi++){
    // //rotation  
    G4RotationMatrix* rotCounter = new G4RotationMatrix;
    rotCounter->rotateZ(isi*60.*deg);
    //    rotCounter->rotateY(isi*60.*deg);
    // //position
    G4ThreeVector xyzCounter(0.,0.,tgt_center_z+50.*mm);
    // //transform
    G4Transform3D posCounter(*rotCounter,xyzCounter);
    
    G4PVPlacement* SiTracker = new G4PVPlacement(posCounter,"SiTracker",SiLV,physiExpHall,FALSE,isi);
    numdet++;
  }

  //region
  // static int init=1;
  // if(init){
  //   fArrayRegion = new G4Region("Array");
  //   SiLV->SetRegion(fArrayRegion);
  //   fArrayRegion->AddRootLogicalVolume(SiLV);
  //   init = 0;
  // }


  //Register to Sensitive Detector

  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  static int init2=1;
  if(init2){
    fArraySD = new SiSD("SiSD",
			"SiCollection",
			numdet);
    SDman->AddNewDetector(fArraySD);

    SiLV->SetSensitiveDetector(fArraySD);
    init2 = 0;
  }


  // Chambers? ======================================================
  {
    // G4Material *Al = materialMgr->GetMaterial("Al");
    // //  G4Box *aBox = new G4Box("BoxA", 10.0*cm, 5.0*cm, 1.0*cm);
    // G4Tubs* singleSolid= new G4Tubs("SINGLE", 
    // 				    4*cm, 4.3*cm, 72.*cm,
    // 				    // 7.92*cm, 15*cm, 72.*cm,
    // 				    0.*degree, 180.*degree );
    // G4LogicalVolume *fLogicBoxA = new G4LogicalVolume(singleSolid,
    // 						      Al,
    // 						      "logicBoxA",
    // 						      0,0,0);
    // G4RotationMatrix* rotCounter = new G4RotationMatrix;
    // rotCounter->rotateY(0.*deg);
    // G4ThreeVector xyzCounter(0.,0.,-4491.*mm);
    // G4Transform3D posCounter(*rotCounter,xyzCounter);
    // G4PVPlacement* target_chamber = new G4PVPlacement(posCounter,"target_chamber",fLogicBoxA,physiExpHall, FALSE, 0);
    
  }
  //minos
  // {
  //   G4Material *P10 = materialMgr->GetMaterial("P10");
  //   //  G4Box *aBox = new G4Box("BoxA", 10.0*cm, 5.0*cm, 1.0*cm);
  //   G4Tubs* singleSolid= new G4Tubs("SINGLE", 
  // 				  4.5*cm, 8.5*cm, 30.*cm,
  // 				  // 7.92*cm, 15*cm, 72.*cm,
  // 				  0.*degree, 180.*degree );
  // G4LogicalVolume *fLogicBoxB = new G4LogicalVolume(singleSolid,
  //                                                       P10,
  //                                                       "logicBoxB",
  //                                                       0,0,0);
  // G4RotationMatrix* rotCounter = new G4RotationMatrix;
  // rotCounter->rotateY(0.*deg);
  // G4ThreeVector xyzCounter(0.,0.,-4600.*mm);
  // G4Transform3D posCounter(*rotCounter,xyzCounter);
  // G4PVPlacement* minos_chamber = new G4PVPlacement(posCounter,"minos_chamber",fLogicBoxB,physiExpHall, FALSE, 0);
  // G4SDManager *SDman = G4SDManager::GetSDMpointer();
  // static int init2=1;

  // fArraySD = new TrackerSD("arraySD",
  // 			       "arrayCollection",
  // 			       1);
  // SDman->AddNewDetector(fLogicBoxB);
  // fLogicBoxB->SetSensitiveDetector(fArraySD);

  // G4SDManager *SDman = G4SDManager::GetSDMpointer();
  // static int init2=1;
  // if(init2){
  //   fArraySD = new CalorimeterSD("arraySD",
  // 				 "arrayCollection",
  // 				 numdet);
  //   SDman->AddNewDetector(fArraySD);

  //   fLogicCrystal1->SetSensitiveDetector(fArraySD);
  //   fLogicCrystal2->SetSensitiveDetector(fArraySD);
  //   fLogicCrystal3->SetSensitiveDetector(fArraySD);
  //   fLogicCrystal4->SetSensitiveDetector(fArraySD);
  //   fLogicCrystal5->SetSensitiveDetector(fArraySD);
  //   fLogicCrystal6->SetSensitiveDetector(fArraySD);
  //   init2 = 0;

// }
  // taget holder ==================================================
  
  // G4Material *pla = materialMgr->GetMaterial("Scinti");

  // G4Tubs *holder_us_body = new G4Tubs("holder_us_body", 
  // 				  4.0*cm, 7.43*cm, 0.15*cm,
  // 				  0.*degree, 360.*degree );
  // G4Tubs *holder_us_body_sub = new G4Tubs("holder_us_body_sub", 
  // 				  5.0*cm, 7.5*cm, 0.151*cm,
  // 				  -40.*degree, 260.*degree );
  // G4SubtractionSolid *TgtHolder_us=new G4SubtractionSolid("TgtHolder_us", holder_us_body, holder_us_body_sub);
  // G4LogicalVolume *holder_us_LV = new G4LogicalVolume(TgtHolder_us,
  //                                                  pla,
  //                                                  "holder_us_LV",
  //                                                  0,0,0);

  // G4Tubs *holder_ds_body = new G4Tubs("holder_ds_body", 
  // 				  4.0*cm, 7.43*cm, 0.1*cm,
  // 				  0.*degree, 360.*degree );
  // G4Tubs *holder_ds_body_sub = new G4Tubs("holder_ds_body_sub", 
  // 				  5.0*cm, 7.5*cm, 0.101*cm,
  // 				  -40.*degree, 260.*degree );
  // G4SubtractionSolid *TgtHolder_ds=new G4SubtractionSolid("TgtHolder_ds", holder_ds_body, holder_ds_body_sub);
  // G4LogicalVolume *holder_ds_LV = new G4LogicalVolume(TgtHolder_ds,
  //                                                  pla,
  //                                                  "holder_ds_LV",
  //                                                  0,0,0);

  // // holder leg----------------------------
  // G4Tubs *holder_leg = new G4Tubs("holder_leg", 
  // 				  5.43*cm, 7.93*cm, 5.*cm,
  // 				  -51.*degree, 2.*degree );
  // G4LogicalVolume *holder_leg_LV = new G4LogicalVolume(holder_leg,
  //                                                  pla,
  //                                                  "holder_leg_LV",
  //                                                  0,0,0);
  // G4Tubs *holder_leg2 = new G4Tubs("holder_leg2", 
  // 				  5.43*cm, 7.93*cm, 5.*cm,
  // 				  -131.*degree, 2.*degree );
  // G4LogicalVolume *holder_leg2_LV = new G4LogicalVolume(holder_leg2,
  //                                                  pla,
  //                                                  "holder_leg2_LV",
  //                                                  0,0,0);


  // // for test ------------------------------------------------------
  // G4Trap *test_box_sol=new G4Trap("test_box_sol",2*cm,0*deg,30*deg,1*cm,1*cm,1*cm,0*deg,1*cm,1*cm,1*cm,0*deg);
  // //  G4Box *test_box_sol=new G4Box("test_box_sol",0.2*cm,50.0*cm,0.2*cm);
  // G4LogicalVolume *test_box_LV = new G4LogicalVolume(test_box_sol,
  //                                                  pla,
  //                                                  "test_box_LV",
  //                                                  0,0,0);
  
  


  //  G4double tgt_center_z= 26.2*mm;

  // //for source run ===============================================

  // G4RotationMatrix* rotCounter_tgtholder = new G4RotationMatrix;
  // rotCounter_tgtholder->rotateY(0.*deg);
  // G4double tgtholder_center_z=tgt_center_z-0.25*cm;//--------------------

  // G4ThreeVector tgtholder_us_xyz(0.,0.,tgtholder_center_z-0.15*cm);
  // G4Transform3D tgtholder_us_pos(*rotCounter_tgtholder,tgtholder_us_xyz);
  // G4ThreeVector tgtholder_ds_xyz(0.,0.,tgtholder_center_z+0.1*cm);
  // G4Transform3D tgtholder_ds_pos(*rotCounter_tgtholder,tgtholder_ds_xyz);

  // G4PVPlacement* TgtHolder_us_BodyPVP = new G4PVPlacement(tgtholder_us_pos,"TgtHolder_us_BodyPVP",holder_us_LV ,physiExpHall, FALSE, 0);
  // G4PVPlacement* TgtHolder_ds_BodyPVP = new G4PVPlacement(tgtholder_ds_pos,"TgtHolder_ds_BodyPVP",holder_ds_LV ,physiExpHall, FALSE, 0);


  // G4RotationMatrix* tgtholder_leg_rot = new G4RotationMatrix;
  // tgtholder_leg_rot->rotateZ(0.*deg);
  // G4ThreeVector tgtholder_leg_xyz(0.,0.,tgtholder_center_z);
  // G4Transform3D tgtholder_leg_pos(*tgtholder_leg_rot,tgtholder_leg_xyz);
  // G4PVPlacement* TgtHolder_leg_PVP = new G4PVPlacement(tgtholder_leg_pos,"TgtHolder_leg_PVP",holder_leg_LV ,physiExpHall, FALSE, 0);
  // G4PVPlacement* TgtHolder_leg2_PVP = new G4PVPlacement(tgtholder_leg_pos,"TgtHolder_leg2_PVP",holder_leg2_LV ,physiExpHall, FALSE, 0);


  //  G4PVPlacement* test_box_PVP = new G4PVPlacement(tgtholder_us_pos,"test_box_PVP",test_box_LV ,physiExpHall, FALSE, 0);
  

  //for beam run ===============================================


  // Materials ======================================================
  // target
  G4double AH= 1.0*g/mole;
  G4double AC= 12.0*g/mole;
  G4Element* elH=new G4Element("Hidrogen","N",1,AH);
  G4Element* elC=new G4Element("Carbon","C",6,AC);

  G4Material* polyeth=new G4Material("polyeth",0.94*g/cm3, 2);
  polyeth->AddElement(elH,2);
  polyeth->AddElement(elC,1);

  G4Box *targetbox = new G4Box("targetbox", 4/2*cm, 4/2*cm, 0.1/2*mm);
  G4LogicalVolume *logictarget = new G4LogicalVolume(targetbox, polyeth, "logictarget");

  G4RotationMatrix* rottgt = new G4RotationMatrix;
  rottgt->rotateY(0.*deg);
//  rottgt->rotateY(-45.*deg);


  G4ThreeVector position(0,0,tgt_center_z);
  G4Transform3D targetpos1(*rottgt,position);
  G4PVPlacement *phystarget = new G4PVPlacement(targetpos1, 
						"phystarget1",
						logictarget,
						physiExpHall,
						FALSE,
						1000);

  // for(G4int i=0;i<11;i++){
  //   G4ThreeVector position(0,0,tgt_center_z+i*5.*cm-10*cm);
  //   G4Transform3D targetpos1(*rottgt,position);
  //   G4PVPlacement *phystarget = new G4PVPlacement(targetpos1, 
  // 						  Form("phystarget%d",i),
  // 						  logictarget,
  // 						  physiExpHall,
  // 						  FALSE,
  // 						  i*1000);
  // }


  /*


  G4Material *Pb = materialMgr -> GetMaterial("Pb");
  G4LogicalVolume *logicTMP = new G4LogicalVolume(tmpbox, Pb, "logictmp");
  
  G4PVPlacement *physiTMP = new G4PVPlacement(0, 
					      G4ThreeVector(),
					      "TMP_phys",
					      logicTMP,
					      physiExpHall,
					      false,
					      0);
  
  G4VisAttributes *tmpVisAtt = new G4VisAttributes(G4Color(1,0,1));
  logicTMP -> SetVisAttributes(tmpVisAtt);
  
  */
  return physiExpHall;

}

// UpdateGeometry ===========================================================
void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}


