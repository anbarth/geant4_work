#ifndef RunAction_hh
#define RunAction_hh

#include "G4UserRunAction.hh"

class DetectorConstruction;
class Run;

class RunAction : public G4UserRunAction
{
public:
    explicit RunAction(const DetectorConstruction* det);
    ~RunAction() override = default;

    G4Run* GenerateRun() override;
    void EndOfRunAction(const G4Run* run) override;

private:
    // Writes spectrum.csv: electrons crossing the spectrum plane, per decay
    void WriteSpectrum(const Run* run, G4int nEvents) const;

    const DetectorConstruction* fDetector;

    // Used only to convert dose per decay into a dose rate
    static constexpr double kActivity_Bq = 1.0e6;
};

#endif
