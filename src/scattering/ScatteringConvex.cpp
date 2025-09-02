#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
                                   bool isOpticalPath, int nActs)
    : Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

bool ScatteringConvex::ScatterLight(double beta, double gamma,
                                    std::vector<Beam> &outBeams)
{
    // m_particle->Rotate(beta, gamma, 0);

    m_incidentEnergy = 0;
    m_treeSize = 0;

    /// first extermal beam
    for (int facetID = 0; facetID < m_particle->nFacets; ++facetID)
    {
        const Point3f &inNormal = m_facets[facetID].in_normal;
        splitting.ComputeCosA(m_incidentDir, inNormal);

        if (!splitting.IsIncident()) /// beam is not incident to this facet
        {
            continue;
        }

        Beam inBeam, outBeam;
        SplitLightToBeams(facetID, inBeam, outBeam);

        auto newId = RecomputeTrackId(0, facetID);

        outBeam.id = newId;
        outBeam.lastFacetId = facetID;
        outBeam.nActs = 0;
        outBeams.push_back(outBeam);

        inBeam.id = newId;
        PushBeamToTree(inBeam, facetID, 0, Location::In);

#ifdef _CHECK_ENERGY_BALANCE
        ComputeFacetEnergy(facetID, outBeam);
#endif
    }

    TraceInternalBeams(outBeams);
    return true;
}

bool ScatteringConvex::ScatterLight(double, double, const std::vector<std::vector<int>> &/*tracks*/, std::vector<Beam> &)
{
    return true;
}

void ScatteringConvex::TraceInternalBeams(std::vector<Beam> &outBeams)
{
    while (m_treeSize != 0)
    {
        Beam beam = m_beamTree[--m_treeSize];

        if (IsTerminalAct(beam))
        {
            continue;
        }

        for (int id = 0; id < m_particle->nFacets; ++id)
        {
            if (id == beam.lastFacetId)
            {
                continue;
            }

            Beam inBeam;
            bool isIncident = SplitSecondaryBeams(beam, id, inBeam, outBeams);

            if (!isIncident)
            {
                continue;
            }

            inBeam.id = RecomputeTrackId(beam.id, id);
            inBeam.locations = beam.locations;
            PushBeamToTree(inBeam, id, beam.nActs+1, Location::In);
        }
    }
}

bool ScatteringConvex::SplitSecondaryBeams(Beam &incidentBeam, int facetID,
                                           Beam &inBeam, std::vector<Beam> &outBeams)
{
    Beam outBeam;
    const Point3f &incidentDir = incidentBeam.direction;

    // ext. normal uses in this calculating
    const Point3f &normal = m_facets[facetID].ex_normal;
    splitting.ComputeCosA(normal, incidentDir);

    if (!splitting.IsIncident())
    {
        return false;
    }

    Intersect(facetID, incidentBeam, outBeam);

    if (outBeam.nVertices < MIN_VERTEX_NUM)
    {
        return false;
    }

#ifdef _DEBUG // DEB
    auto id = RecomputeTrackId(incidentBeam.id, facetID);
    if (id == 423)
        int fff = 0;
#endif
    inBeam = outBeam;
    splitting.ComputeSplittingParams(incidentBeam.direction, normal);
//    ComputePolarisationParams(incidentBeam.direction, normal, incidentBeam);

    if (!splitting.IsNormalIncidence())
    {	// regular incidence
        Beam incBeam = incidentBeam;

        ComputePolarisationParams(incidentBeam.direction, normal, incBeam);

        if (!splitting.IsCompleteReflection())
        {
            outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);

            splitting.ComputeRegularBeamsParams(normal, incBeam, inBeam, outBeam);
            outBeam.nActs = incidentBeam.nActs + 1;
            outBeam.opticalPath += splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
            outBeam.lastFacetId = facetID;
            outBeams.push_back(outBeam);
        }
        else // complete internal reflection incidence
        {
            splitting.ComputeCRBeamParams(normal, incBeam, inBeam);
        }
    }
    else
    {	// normal incidence
        splitting.ComputeNormalBeamParams(incidentBeam, inBeam, outBeam);

        outBeam.nActs = incidentBeam.nActs + 1;
        outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);
        double path = splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
        outBeam.opticalPath += path;
        outBeam.lastFacetId = facetID;
        outBeams.push_back(outBeam);
    }

    return true;
}
