#pragma once

#include "Handler.h"

class HandlerPO : public Handler
{
public:
    HandlerPO(Particle *particle, Light *incidentLight, int nTheta,
              double wavelength);

    void HandleBeams(std::vector<Beam> &beams, double sinZenith) override;
    void WriteMatricesToFile(std::string &destName, double nrg) override;
    void WriteTotalMatricesToFile(const std::string &destName) override;
    // double ComputeTotalScatteringEnergy() override;

    void SetScatteringSphere(const ScatteringRange &grid) override;

    void SetBackScatteringConus(double radAngle);

    matrix *m_Lp;
    matrix *m_Ln;
    Arr2D M;				// Mueller matrices

protected:
    virtual void AddToMueller();

    void ComputeOpticalLengths(const Beam &beam, BeamInfo &info);

    virtual void RotateJones(const Beam &beam, const BeamInfo &info,
                     const Vector3d &vf, const Vector3d &direction,
                     matrixC &matrix) const;
    void CleanJ();
    matrixC ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
                           const Vector3d &direction);

    matrixC ApplyDiffraction(const Beam &beam, const BeamInfo &info,
                         const Vector3d &direction, const Vector3d &vf);

    BeamInfo ComputeBeamInfo(Beam &beam);



protected:
    std::vector<Arr2D> m_groupMatrices;	//
    std::vector<Arr2DC> m_diffractedMatrices;	// Jones matrices
    bool isNanOccured = false;
    bool isNan = false;
    bool isBackScatteringConusEnabled = false;
    double backScatteringConus = 180;

    // Handler interface
public:
    void SetTracks(Tracks *tracks) override;
private:
    void WriteGroupMatrices(Arr2D &matrices, const std::string &name);
};

