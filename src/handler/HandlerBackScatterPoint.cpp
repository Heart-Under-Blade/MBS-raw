#include "HandlerBackScatterPoint.h"

//#define BEAM_DIR_LIM 0.9396 // cos(20)
//#define BEAM_DIR_LIM 0.98480775301 // cos(10)
#define BEAM_DIR_LIM 0.76604444311 // cos(40)

HandlerBackScatterPoint::HandlerBackScatterPoint(Particle *particle,
                                                 Light *incidentLight,
                                                 int nTheta,
                                                 double wavelength)
    : HandlerPO(particle, incidentLight, nTheta, wavelength)
{
}

void HandlerBackScatterPoint::HandleBeams(std::vector<Beam> &beams, double sinZenith)
{
    Point3d backDirection(0, 0, 1);
    Point3d vf = -m_incidentLight->polarizationBasis;

#ifdef _DEBUG // DEB
    int cc = 0;
#endif
    for (Beam &beam : beams)
    {
#ifdef _DEBUG // DEB
        cc++;
//		double m11 = real(diffractedMatrix[0][0]);
        std::vector<int> tr;
        Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
        if (tr[0] == 0 && tr[1] == 1 && tr[2] == 7 && tr[3] == 0)
        {
            int ddfdfdf = 0;
        }
#endif
        if (isBackScatteringConusEnabled && beam.direction.cz < backScatteringConus)
        {
            continue;
        }

        int groupId = m_tracks->FindGroupByTrackId(beam.id);

//        if (groupId < 0 && m_tracks->shouldComputeTracksOnly)
//        {
//            continue;
//        }

#ifdef _DEBUG // DEB
        if (cc == 93)
            int ddddddd = 0;
#endif
        beam.polarizationBasis = beam.RotateSpherical(-m_incidentLight->direction,
                                                      m_incidentLight->polarizationBasis);

        BeamInfo info = ComputeBeamInfo(beam);

        if (m_isBadBeam)
        {
            continue;
        }

        if (beam.lastFacetId != __INT_MAX__)
        {
            ApplyAbsorption(beam);
        }

        matrixC diffractedMatrix = ApplyDiffraction(beam, info, backDirection, vf);
        /* DEB */ {
#ifdef _DEBUG
            complex ddd = diffractedMatrix[0][0];
            int fff = 0;
#endif
        }
        // correction
        Matrix2x2c jonesCor = diffractedMatrix;
        jonesCor.m12 -= jonesCor.m21;
        jonesCor.m12 /= 2;
        jonesCor.m21 = -jonesCor.m12;

//        if (groupId < 0 && !m_tracks->shouldComputeTracksOnly) // tracks not from list
//        {
//            originContrib->AddToMueller(diffractedMatrix);
//            correctedContrib->AddToMueller(jonesCor);
//        }
//        else // tracks from list (coh)
//        {
            originContrib->AddToGroup(diffractedMatrix, groupId);
            correctedContrib->AddToGroup(jonesCor, groupId);
//        }
    }

    originContrib->SumGroupTotal();
    correctedContrib->SumGroupTotal();
}

void HandlerBackScatterPoint::SetTracks(Tracks *tracks)
{
    Handler::SetTracks(tracks);
    originContrib = new PointContribution(tracks->size()+1, m_normIndex);
    correctedContrib = new PointContribution(tracks->size()+1, m_normIndex);
}

void HandlerBackScatterPoint::OutputContribution(ScatteringFiles &files,
                                                 double angle, double energy,
                                                 bool isOutputGroups,
                                                 std::string prefix)
{
    PointContribution *contrib = (prefix == "") ? originContrib : correctedContrib;
    contrib->SumTotal();

    energy *= m_normIndex;

    auto m = contrib->GetTotal();

    std::ofstream *all = files.GetMainFile(prefix + "all");
    *(all) << angle << ' ' << energy << ' ';
    *(all) << m << std::endl;
//cout << endl << endl << contrib->GetRest()(0,0) << endl << endl ;
    if (isOutputGroups)
    {
        for (size_t gr = 0; gr < m_tracks->size(); ++gr)
        {
            if (m_tracks->at(gr).size > 0)
            {
                std::ofstream &file = *(files.GetGroupFile(gr));
                file << angle << ' ' << energy << ' ';
                file << contrib->GetGroupMueller(gr) << std::endl;
            }
        }
    }

    if (!m_tracks->shouldComputeTracksOnly)
    {
        std::ofstream &other = *(files.GetMainFile(prefix + "other"));
        other << angle << ' ' << energy << ' ';
        other << contrib->GetRest() << std::endl;

        std::ofstream &diff = *(files.GetMainFile(prefix + "difference"));
        diff << angle << ' ' << energy << ' ';
        diff << contrib->GetGroupTotal() << std::endl;
    }

    contrib->Reset();
}


void HandlerBackScatterPoint::RotateJones(
        const Beam &beam, const BeamInfo &info, const Vector3d &vf,
        const Vector3d &direction, matrixC &matrix) const
{
    auto &dir = direction;

    Vector3d vt = CrossProductD(vf, dir);
    vt = vt/LengthD(vt);

    Vector3f NT = CrossProduct(info.normal, info.beamBasis);
    Vector3f NE = CrossProduct(info.normal, beam.polarizationBasis);

    Vector3d NTd = Vector3d(NT.cx, NT.cy, NT.cz);
    Vector3d NEd = Vector3d(NE.cx, NE.cy, NE.cz);


    matrix[0][0] = -DotProductD(NTd, vf);
    matrix[0][1] = -DotProductD(NEd, vf);
    matrix[1][0] = DotProductD(NTd, vt);
    matrix[1][1] = DotProductD(NEd, vt);
}
