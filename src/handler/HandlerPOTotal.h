#pragma once

#include "HandlerPO.h"

class HandlerPOTotal : public HandlerPO
{
public:
    HandlerPOTotal(Particle *particle, Light *incidentLight, int nTheta,
                   double wavelength);
    ~HandlerPOTotal() { delete betaMueller; }
    void WriteMatricesToFile(std::string &destName, double nrg) override;
    void AddToMueller() override;

    void OutputContribution(double angle, double energy);

    matrix *betaMueller;
};
