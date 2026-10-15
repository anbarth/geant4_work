#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1120
#include "G4HadronicParameters.hh"
#endif

#include "PMDetectorConstruction.hh"
#include "PMActionInitialization.hh"

int main(int argc, char** argv)
{
    // Geant4 >= 11.2 ignores radioactive decays that happen later than a
    // time threshold (default 1 year). Sr-90 has a 28.8 y half-life, so the
    // threshold must be raised or the source would never decay.
#if G4VERSION_NUMBER >= 1120
    G4HadronicParameters::Instance()->SetTimeThresholdForRadioactiveDecay(1.0e+60 * year);
#endif

    // Interactive session only if no macro was given
    G4UIExecutive* ui = nullptr;
    if (argc == 1) ui = new G4UIExecutive(argc, argv);

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    runManager->SetUserInitialization(new PMDetectorConstruction());

    // Physics: precise EM (needed for low-energy betas in thin layers),
    // particle decays, and radioactive decay (Sr90 -> Y90 -> Zr90)
    auto* physics = new G4VModularPhysicsList();
    physics->RegisterPhysics(new G4EmStandardPhysics_option4());
    physics->RegisterPhysics(new G4DecayPhysics());
    physics->RegisterPhysics(new G4RadioactiveDecayPhysics());
    physics->SetDefaultCutValue(0.01 * mm);
    runManager->SetUserInitialization(physics);

    runManager->SetUserInitialization(new PMActionInitialization());

    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    G4UImanager* uiManager = G4UImanager::GetUIpointer();

    if (ui) {
        uiManager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    } else {
        G4String command = "/control/execute ";
        uiManager->ApplyCommand(command + argv[1]);
    }

    delete visManager;
    delete runManager;
    return 0;
}
