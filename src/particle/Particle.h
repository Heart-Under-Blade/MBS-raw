#pragma once

#include <math.h>

#include "compl.hpp"
#include "geometry_lib.h"
#include "Facet.h"
#include <vector>

#define ROT_MTR_RANK 3

class Angle
{
public:
    double alpha;
    double beta;
    double gamma;
};

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle
{
public:
    Particle();

    void SetFromFile(const std::string &filename);

    void Rotate(double beta, double gamma, double alpha);
    void Move(float dx, float dy, float dz);
    void Fix();
    void Resize(double size);

    Facet CreateBase(int nVertices, double radius, double zCoordinate);
    void CreateSideFacets(const Facet &top, const Facet &bottom, int &iFacet);

    void Concate(const std::vector<Particle> &parts);

    /**
     * @brief LongRadius
     * @return The distance from beginning of the center of coordinate system
     * to the farthest point of particle.
     */
    double LongRadius() const;

    /**
     * @brief A sum of areas of each facet of the particle
     * @return value of area
     */
    double Area();

    double Volume();
    Point3f Center() const;

    double MaximalDimention() const;
    double MaximalDimentionPart() const;

    const complex &GetRefractiveIndex() const;
    void SetRefractiveIndex(const complex &value);

    const Symmetry &GetSymmetry() const;
    virtual void GetParticalFacetIdRangeByFacetId(
            int /*id*/, int &/*begin*/, int &/*end*/) const {}

    bool IsConcave() const;

    void Output(std::string name);
    void Scale(double ratio);
    // std::vector<unsigned> FindFlippedFacets() const;
    // void DisconnectParts(std::vector<Particle> &parts) const; ///< Get Particle's the aggregate consist of.

public:
    Facet facets[MAX_FACET_NUM];	///< all facets of particle
    int nFacets;					///< number of facets
    bool isAggregated = false;
    bool isConcave;
    int nFacetsInPart = 0;

    Angle rotAngle;

    Symmetry m_symmetry;		///< angle of particle symmetry
protected:
    Facet defaultFacets[MAX_FACET_NUM];


    complex m_refractiveIndex;	///< complex value of refractive index of the particle

protected:
    void Init(int facetCount, const complex &refrIndex);

    void SetDefaultNormals();
    void SetDefaultCenters();
    void Reset();
    void SetSymmetry(double beta, double gamma, double alpha = 0);

    virtual void SetFacetParams() {}

private:
    void SetDParams();
    void RotateNormals();
    void RotatePoint(const Point3f &point, Point3f &result);
    void RotateCenters();
    void SetRotateMatrix(double beta, double gamma, double alpha);
    void RemoveWalls();
    std::vector<Facet> ConnectBases(const Facet &top, const Facet &bottom);
    void FindFlippedFacets(std::vector<unsigned> &facetIndices) const;

private:
    double m_rotMatrix[ROT_MTR_RANK][ROT_MTR_RANK];	///< rotation matrix for vertices
    void ReadSymmetry(const int bufSize, char *trash, char *buff,
                      std::ifstream pfile, char *ptr);
};

