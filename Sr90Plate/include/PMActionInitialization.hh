#ifndef PMActionInitialization_hh
#define PMActionInitialization_hh

#include "G4VUserActionInitialization.hh"

class PMActionInitialization : public G4VUserActionInitialization
{
public:
    PMActionInitialization() = default;
    ~PMActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;
};

#endif
