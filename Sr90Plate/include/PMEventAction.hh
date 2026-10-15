#ifndef PMEventAction_hh
#define PMEventAction_hh

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

class PMEventAction : public G4UserEventAction
{
public:
    PMEventAction();
    ~PMEventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void AddEdep(G4int well, G4double edep) { fEdep[well] += edep; }

private:
    std::vector<G4double> fEdep;
};

#endif
