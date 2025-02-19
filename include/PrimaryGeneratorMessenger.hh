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
// $Id: PrimaryGeneratorMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PRIMARYGENERATORMESSENGER_HH
#define PRIMARYGENERATORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
class PrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger : public G4UImessenger
{
public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction *);
    virtual ~PrimaryGeneratorMessenger();

public:
    void SetNewValue(G4UIcommand *, G4String);

private:
    PrimaryGeneratorAction *fAction;

    G4UIcmdWithAnInteger *fNumberOfDecayingLaBrDetectorsCmd;
    G4UIcmdWithADoubleAndUnit *fEfficiencyEnergyCmd;
    G4UIcmdWith3Vector *fEfficiencyDirectionCmd;
    G4UIcmdWith3VectorAndUnit *fEfficiencyPositionCmd;
    G4UIcmdWithAString *fEfficiencyParticleCmd;
    G4UIcmdWith3Vector *fEfficiencyPolarizationCmd;
    G4UIcmdWithADoubleAndUnit *fConeRadiusCmd; // SPICE cone radius
    G4UIcmdWithADoubleAndUnit *fConeZValueCmd; // SPICE cone with values
    G4UIcmdWithADoubleAndUnit *fConeRValueCmd;
    G4UIcmdWithADoubleAndUnit *fConeAngleCmd;
    G4UIcmdWithADoubleAndUnit *fConeMinAngleCmd;
    G4UIcmdWithADoubleAndUnit *fBeamSpotSigmaCmd;
    G4UIcmdWithADoubleAndUnit *fSourceRadiusCmd;
    G4UIcmdWithAnInteger *fBeamDistroCmd;
    G4UIcmdWithAString *fBeamFileCmd;
    G4UIcmdWithADoubleAndUnit *fKentuckyEnergyCmd;
    G4UIcmdWithAString *fKentuckyReactionCmd;
    G4UIcmdWithADoubleAndUnit *fMinimumPhiCmd;
    G4UIcmdWithADoubleAndUnit *fMaximumPhiCmd;
    G4UIcmdWithADoubleAndUnit *fMinimumThetaCmd;
    G4UIcmdWithADoubleAndUnit *fMaximumThetaCmd;
    G4UIcmdWithAnInteger *fVerbosityCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
