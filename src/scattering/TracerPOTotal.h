#pragma once

#include "TracerPO.h"

class TracerPOTotal : public TracerPO
{
public:
	TracerPOTotal(Particle *particle, int nActs, const std::string &resultFileName);
	void TraceRandom(const AngleRange &betaRange,
					 const AngleRange &gammaRange) override;
    void TraceMonteCarlo(const AngleRange &betaRange,
                         const AngleRange &gammaRange, int nOrientations);
};
