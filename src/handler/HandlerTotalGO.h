#pragma once

#include "HandlerGO.h"

class HandlerTotalGO : public HandlerGO
{
public:
    HandlerTotalGO(Particle *particle, Light *incidentLight, int nTheta,
                   float wavelength);

    void HandleBeams(std::vector<Beam> &beams, double sinZenith) override;
    void WriteMatricesToFile(std::string &destName, double nrg) override;

    void SetScatteringSphere(const ScatteringRange &grid) override;
};

