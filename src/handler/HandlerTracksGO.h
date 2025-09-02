#pragma once

#include "HandlerGO.h"

class HandlerTracksGO : public HandlerGO
{
public:
    HandlerTracksGO(Particle *particle, Light *incidentLight, int nTheta,
                    float wavelength);

    void HandleBeams(std::vector<Beam> &beams, double sinZenith) override;
    void WriteMatricesToFile(std::string &destName, double nrg) override;
};
