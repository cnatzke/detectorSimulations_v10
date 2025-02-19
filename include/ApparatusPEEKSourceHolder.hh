//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ApparatusPEEKSourceHolder_h
#define ApparatusPEEKSourceHolder_h 1

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4AssemblyVolume;
class G4Material;

class ApparatusPEEKSourceHolder
{
public:
    ApparatusPEEKSourceHolder();
    ~ApparatusPEEKSourceHolder();

    G4int Build(); // G4SDManager* mySDman);
    G4int Place(G4LogicalVolume *expHallLog, G4int detectorNumber);
    G4double GetR();
    G4double GetTheta(G4int i);
    G4double GetPhi(G4int i);
    G4double GetYaw(G4int i);
    G4double GetPitch(G4int i);
    G4double GetRoll(G4int i);

private:
    // Material Manager
    G4NistManager *nistMan;

    // Logical volumes
    G4LogicalVolume *fCeramicPelletLog;
    G4LogicalVolume *fSourceSphereLog;
    G4LogicalVolume *fSourceSupportLog;

    // Assembly volumes
    G4AssemblyVolume *fAssembly;

    // Surface Check
    G4bool fSurfCheck;

    // Solid Volumes
    G4Sphere *sourceSphere;
    G4Tubs *cyc;
    G4SubtractionSolid *delSphere;

    // Parameters

    G4int fCopyNumber;
    G4int fNumberOfDetectors;
    G4int fNumberOfSegments;

    G4String fSphereMaterial;
    G4String fSupportMaterial;

    G4Material *fPelletMaterial;

    G4double fullPhi;
    G4double fullTheta;
    G4double fPelletInnerRadius;
    G4double fPelletOuterRadius;
    G4double fPelletLength;
    G4double fSphereInnerRadius;
    G4double fSphereOuterRadius;
    G4double fLargeCyclinderLength;
    G4double fLargeCyclinderInnerRadius;
    G4double fLargeCyclinderOuterRadius;
    G4double fMediumCyclinderLength;
    G4double fMediumCyclinderInnerRadius;
    G4double fMediumCyclinderOuterRadius;
    G4double fSmallCyclinderLength;
    G4double fSmallCyclinderInnerRadius;
    G4double fSmallCyclinderOuterRadius;

    G4double fDetectorAngles[8][5];

    G4Tubs *BuildCeramicPellet();
    G4Sphere *BuildSourceSphere();
    G4Tubs *BuildLargeSourceCyclinder();
    G4Tubs *BuildMediumSourceCyclinder();
    G4Tubs *BuildSmallSourceCyclinder();

    G4int BuildCeramicPelletVolume();
    G4int BuildSourceSphereVolume();
    G4int BuildSourceSupportVolume();

    G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);
};

#endif
