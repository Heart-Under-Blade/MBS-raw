#pragma once

#include "HandlerPO.h"
#include "ScatteringFiles.h"

class HandlerBackScatterPoint : public HandlerPO
{
public:
    HandlerBackScatterPoint(Particle *particle, Light *incidentLight, int nTheta,
                            double wavelength);

    void HandleBeams(std::vector<Beam> &beams, double sinZenith) override;
    void SetTracks(Tracks *tracks) override;

    void OutputContribution(ScatteringFiles &files, double angle, double energy,
                            bool isOutputGroups, std::string prefix = "");

    PointContribution *originContrib;
private:
    PointContribution *correctedContrib;


    // HandlerPO interface
protected:
    void RotateJones(const Beam &beam, const BeamInfo &info, const Vector3d &vf, const Vector3d &direction, matrixC &matrix) const override;
};

