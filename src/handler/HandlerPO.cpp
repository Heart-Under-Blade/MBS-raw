#include "HandlerPO.h"

#include "Mueller.hpp"
#include <iostream>

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, int nTheta,
                     double wavelength)
    : Handler(particle, incidentLight, nTheta, wavelength)
{
    m_Lp = new matrix(4, 4);

    (*m_Lp)[0][0] = 1;
    (*m_Lp)[0][1] = 0;
    (*m_Lp)[0][2] = 0;
    (*m_Lp)[0][3] = 0;

    (*m_Lp)[1][0] = 0;
    (*m_Lp)[1][3] = 0;

    (*m_Lp)[2][0] = 0;
    (*m_Lp)[2][3] = 0;

    (*m_Lp)[3][0] = 0;
    (*m_Lp)[3][1] = 0;
    (*m_Lp)[3][2] = 0;
    (*m_Lp)[3][3] = 1;

    m_Ln = new matrix(4, 4);

    (*m_Ln) = (*m_Lp);
}

void HandlerPO::CleanJ()
{
    m_diffractedMatrices.clear();
    Arr2DC tmp(m_sphere.nAzimuth + 1, m_sphere.nZenith + 1, 2, 2);
    tmp.ClearArr();

    for (unsigned q = 0; q < m_tracks->size() + 1; q++)
    {
        m_diffractedMatrices.push_back(tmp);
    }
}

void HandlerPO::WriteMatricesToFile(std::string &destName, double nrg)
{
//    if (!m_tracks->shouldComputeTracksOnly)
    {
        std::ofstream outFile(destName, std::ios::app);

        outFile /*<< std::to_string(m_sphere.radius) << ' '*/
                << std::to_string(m_sphere.nZenith) << ' '
                << std::to_string(m_sphere.nAzimuth+1);

        for (int t = 0; t <= m_sphere.nZenith; ++t)
        {
            double tt = RadToDeg(t*m_sphere.zenithStep);

            for (int p = 0; p <= m_sphere.nAzimuth; ++p)
            {
                double fi = -((double)p)*m_sphere.azinuthStep;
                double degPhi = RadToDeg(-fi);
                outFile << std::endl << tt << " " << degPhi << " " << nrg << " ";

                matrix m = M(p, t);
                outFile << m;
            }
        }
    }
}

void HandlerPO::WriteGroupMatrices(Arr2D &matrices, const std::string &name)
{
    auto &Lp = *m_Lp;
    auto &Ln = *m_Ln;

    std::ofstream outFile(name, std::ios::out);

    outFile /*<< std::to_string(m_sphere.radius) << ' '*/
            << std::to_string(m_sphere.nZenith) << ' '
            << std::to_string(m_sphere.nAzimuth+1);

    matrix sum(4, 4);

    for (int t = m_sphere.nZenith; t >= 0; --t)
    {
        sum.Fill(0.0);
        double tt = RadToDeg(m_sphere.zenithEnd-m_sphere.zenithStart) - RadToDeg(t*m_sphere.zenithStep);

        for (int p = 0; p <= m_sphere.nAzimuth; ++p)
        {
            double fi = -((double)p)*m_sphere.azinuthStep;
            matrix m = matrices(p, t);

            Lp[1][1] = cos(2*fi);
            Lp[1][2] = sin(2*fi);
            Lp[2][1] = -Lp[1][2];
            Lp[2][2] = Lp[1][1];

            Ln[1][2] = -Lp[1][2];
            Ln[2][1] = -Lp[2][1];

            if (t == 0)
            {
                sum += Lp*m*Lp;
            }
            else if (t == m_sphere.nZenith-1)
            {
                sum += Ln*m*Lp; // OPT: вынести Ln в отдельный случай
            }
            else
            {
                sum += m*Lp;
            }
        }

        outFile << std::endl << tt << " ";
        outFile << sum/m_sphere.nAzimuth;
    }
}

void HandlerPO::WriteTotalMatricesToFile(const std::string &destName)
{
    if (m_tracks->shouldOutputGroups)
    {
        for (size_t i = 0; i < m_diffractedMatrices.size(); ++i)
        {
            if ((*m_tracks)[i].size != 0)
            {
                const std::string subname = (*m_tracks)[i].CreateGroupName();
                const std::string &filename = destName + '_' +  subname;
                WriteGroupMatrices(m_groupMatrices[i], filename);
            }
        }

        WriteGroupMatrices(M, destName + "_total.dat");
    }
    std::cout << std::endl << destName;
}

// double HandlerPO::ComputeTotalScatteringEnergy()
// {
//     double D_tot = 0;

//     for (int t = 0; t <= m_sphere.nZenith; ++t)
//     {
//         for (int p = 0; p <= m_sphere.nAzimuth; ++p)
//         {
//             matrix m = M(p, t);
//             D_tot += m[0][0];
//         }
//     }

//     return D_tot;
// }

void HandlerPO::RotateJones(const Beam &beam, const BeamInfo &info,
                            const Vector3d &vf, const Vector3d &direction,
                            matrixC &matrix) const
{
    auto &dir = direction;

    Vector3d vt = CrossProductD(vf, dir);
    vt = vt/LengthD(vt);

    Vector3f NT = CrossProduct(info.normal, info.beamBasis);
    Vector3f NE = CrossProduct(info.normal, beam.polarizationBasis);

    Vector3d NTd = Vector3d(NT.cx, NT.cy, NT.cz);
    Vector3d NPd = Vector3d(NE.cx, NE.cy, NE.cz);

    Point3f DT = CrossProduct(beam.direction, info.beamBasis);
    Point3f DP = CrossProduct(beam.direction, beam.polarizationBasis);

    Point3d DTd = Point3d(DT.cx, DT.cy, DT.cz);
    Point3d DPd = Point3d(DP.cx, DP.cy, DP.cz);

    const Point3d nd = info.normal;

    Point3d cpT = CrossProductD(dir, NTd)
            - CrossProductD(dir, CrossProductD(dir, CrossProductD(nd, DTd)));
    Point3d cpP = CrossProductD(dir, NPd)
            - CrossProductD(dir, CrossProductD(dir, CrossProductD(nd, DPd)));

    matrix[0][0] = DotProductD(cpT, vt)/2.0;
    matrix[0][1] = DotProductD(cpP, vt)/2.0;
    matrix[1][0] = DotProductD(cpT, vf)/2.0;
    matrix[1][1] = DotProductD(cpP, vf)/2.0;
}

matrixC HandlerPO::ApplyDiffraction(const Beam &beam, const BeamInfo &info,
                                    const Vector3d &direction,
                                    const Vector3d &vf)
{
    matrixC fnJones = (beam.lastFacetId != __INT_MAX__) ?
                ComputeFnJones(beam.J, info, direction) :
                beam.J * exp_im(m_waveIndex*beam.opticalPath);

    matrixC jones_rot(2, 2);
    RotateJones(beam, info, vf, direction, jones_rot);

    complex fresnel = DiffractIncline(info, beam, direction);
//	complex fresnel = (m_hasAbsorption && beam.lastFacetId != INT_MAX && beam.nActs > 0)
//			? DiffractInclineAbs(info, beam, direction)
//			: DiffractIncline(info, beam, direction);

    if (isnan(real(fresnel)))
    {
        return matrixC(2, 2);
    }

    matrixC difracted = fresnel*jones_rot*fnJones;
#ifdef _DEBUG // DEB
    complex ddd[4];
    ddd[0] = jones_rot[0][0];
    ddd[1] = jones_rot[0][1];
    ddd[2] = jones_rot[1][0];
    ddd[3] = jones_rot[1][1];

    complex qqq[4];
    qqq[0] = fnJones[0][0];
    qqq[1] = fnJones[0][1];
    qqq[2] = fnJones[1][0];
    qqq[3] = fnJones[1][1];

    complex bbb[4];
    bbb[0] = difracted[0][0];
    bbb[1] = difracted[0][1];
    bbb[2] = difracted[1][0];
    bbb[3] = difracted[1][1];
#endif
    return difracted;
}

matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
                                  const Vector3d &direction)
{
    double dp = DotProductD(direction, info.center);
    double arg = m_waveIndex*(info.projLenght - dp);
    return matrix*exp_im(arg);
}

void HandlerPO::AddToMueller()
{
    for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
    {
        for (int t = 0; t <= m_sphere.nZenith; ++t)
        {
            for (int p = 0; p <= m_sphere.nAzimuth; ++p)
            {
                matrix m = Mueller(m_diffractedMatrices[q](p, t));
                m *= m_normIndex;
                M.insert(p, t, m);
                m_groupMatrices[q].insert(p, t, m);
            }
        }
    }
}

BeamInfo HandlerPO::ComputeBeamInfo(Beam &beam)
{
    BeamInfo info;
    info.normal = beam.Normal();
    info.normald = Point3d(info.normal.cx, info.normal.cy, info.normal.cz);

    info.order = DotProduct(info.normal, beam.direction) > 0;

    if (!info.order)
    {
        //info.normal = -info.normal;
//        info.normald = -info.normald;
    }

    m_isBadBeam = false;

    ComputeCoordinateSystemAxes(info.normald, info.horAxis, info.verAxis);

    info.center = beam.Center();
    info.projectedCenter = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                                  info.normald, info.center);

    if (m_hasAbsorption && beam.lastFacetId != __INT_MAX__ && beam.nActs > 1)
    {
        ComputeOpticalLengths(beam, info);
        ComputeLengthIndices(beam, info);
    }

    info.area = beam.Area();

    info.projLenght = beam.opticalPath + DotProductD(info.center, beam.direction);
//#ifdef _DEBUG // DEB
//	std::vector<int> tr;
//	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

//	double sss = m_scattering->ComputeInternalOpticalPath(
//				beam, beam.Center(), tr);
//#endif
    info.beamBasis = CrossProduct(beam.polarizationBasis, beam.direction);
    info.beamBasis = info.beamBasis/Length(info.beamBasis); // basis of beam

    return info;
}

void HandlerPO::SetBackScatteringConus(double radAngle)
{
    isBackScatteringConusEnabled = true;
    backScatteringConus = cos(radAngle);
}

void HandlerPO::SetScatteringSphere(const ScatteringRange &grid)
{
    m_sphere = grid;
    M = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);

    m_sphere.ComputeSphereDirections(*m_incidentLight);
}

void HandlerPO::ComputeOpticalLengths(const Beam &beam, BeamInfo &info)
{
    std::vector<int> tr;
    Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

    for (int i = 0; i < 3; ++i)
    {
        info.opticalLengths[i] = m_scattering->ComputeInternalOpticalPath(
                    beam, beam.arr[i], tr);
    }

//#ifdef _DEBUG // DEB
//	if (info.opticalLengths[0] > 200)
//		int fff = 0;
//#endif

//	double path = m_scattering->ComputeInternalOpticalPath(beam, beam.Center(), tr);

//	if (path > DBL_EPSILON)
//	{
//		double abs = exp(m_cAbs*path);
//		beam.J *= abs;
//	}

    //	ExtropolateOpticalLenght(beam, tr);
}

void HandlerPO::HandleBeams(std::vector<Beam> &beams, double sinZenith)
{
    m_sinZenith = sinZenith;
    double sum = 0;
    int groupId;
    CleanJ();

    for (Beam &beam : beams)
    {
        if (isBackScatteringConusEnabled && beam.direction.cz < backScatteringConus)
        {
            continue;
        }

        groupId = 0;
#ifdef _DEBUG // DEB
        // if (beam.id == 423)
        //     int dsdsdsd = 0;
//		std::vector<int> tr;
//		Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
//		if (tr.size() == 2 && tr[0] == 2 && tr[1] == 4)
//			int fff = 0;
#endif
//		if (m_tracks->shouldComputeTracksOnly)
//		{
//			groupId = m_tracks->FindGroupByTrackId(beam.id);

            if (groupId < 0)
            {
//				groupId = 0;
//				continue;
            }
//		}

        beam.polarizationBasis = beam.RotateSpherical(
                    -m_incidentLight->direction,
                    m_incidentLight->polarizationBasis);

        BeamInfo info = ComputeBeamInfo(beam);

        if (m_isBadBeam)
        {
            continue;
        }

        if (m_hasAbsorption && beam.lastFacetId != __INT_MAX__ && beam.lastFacetId != -1)
        {
            ApplyAbsorption(beam);
        }

        if (beam.lastFacetId != __INT_MAX__)
        {
            matrix m_ = Mueller(beam.J);
            m_outputEnergy += BeamCrossSection(beam)*m_[0][0]*m_sinZenith;
        }

#ifdef _DEBUG // DEB
        sum += beam.Area();
        for (int i = 0; i <= m_sphere.nAzimuth; ++i)
        {
            for (int j = 0; j <= m_sphere.nZenith; ++j)
            {
#else

        for (int i = 0; i <= m_sphere.nAzimuth; ++i)
        {
            for (int j = 0; j <= m_sphere.nZenith; ++j)
            {
#endif
                Point3d &dir = m_sphere.directions[i][j];
                Point3d &vf = (j == 0) ? m_sphere.vf.back() : m_sphere.vf[i];
                matrixC diffractedMatrix = ApplyDiffraction(beam, info, dir, vf);

//                if (groupId < 0)
                if (!isCoh)
                {
                    matrix m = Mueller(diffractedMatrix);
#ifdef _DEBUG // DEB
                    // double &ddd = m[0][0];
#endif
                    m *= m_sinZenith;
                    // m /= 2; // TODO: костыль
                    M.insert(i, j, m);
                }
                else
                {
                    // matrix m = Mueller(diffractedMatrix);
                    // m_outputEnergy += /*xs**/m[0][0];
                    m_diffractedMatrices[groupId].insert(i, j, diffractedMatrix);

#ifdef _DEBUG // DEB
                    if (i == m_sphere.nAzimuth && j == m_sphere.nZenith)
                        int fff = 0;
                    complex ddd = diffractedMatrix[0][0];
//                    std::cout << beam.opticalPath << " "<< real(ddd) << " " << imag(ddd) << " " << beam.nActs << std::endl;
                    int fff = 0;
#endif
                }
            }
        }
#ifdef _DEBUG // DEB
        complex ddd = m_diffractedMatrices[0](0,0)[0][0];
        int fff = 0;
#endif
    }

    if (isCoh)
    {
        AddToMueller();
    }

#ifdef _DEBUG // DEB
//    double mm = M(28, 20)[0][0];
    int fff = 0;
#endif

//    if (m_tracks->shouldComputeTracksOnly)
//    {
//        AddToMueller();
//    }
}

void HandlerPO::SetTracks(Tracks *tracks)
{
    m_tracks = tracks;

    for (int i = 0; i < m_tracks->size() + 1; ++i)
    {
        Arr2D tmp = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);
        m_groupMatrices.push_back(tmp);
    }
}
