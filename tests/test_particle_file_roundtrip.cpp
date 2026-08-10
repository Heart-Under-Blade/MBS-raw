#include "particle/Bullet.h"
#include "particle/BulletRosette.h"
#include "particle/CertainAggregate.h"
#include "particle/ConcaveHexagonal.h"
#include "particle/Droxtal.h"
#include "particle/Hexagonal.h"
#include "particle/HexagonalAggregate.h"
#include "global.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

namespace
{
struct Differences
{
    int topology = 0;
    int coordinates = 0;
    int visibility = 0;
};

Differences CountDifferences(const Particle &expected, const Particle &actual)
{
    Differences differences;
    if (expected.nFacets != actual.nFacets)
    {
        differences.topology = 1000000;
        return differences;
    }
    for (int i = 0; i < expected.nFacets; ++i)
    {
        differences.visibility += expected.facets[i].isVisibleIn != actual.facets[i].isVisibleIn;
        differences.visibility += expected.facets[i].isVisibleOut != actual.facets[i].isVisibleOut;
        differences.topology += expected.facets[i].nVertices != actual.facets[i].nVertices;
        for (int j = 0; j < expected.facets[i].nVertices; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                differences.coordinates += expected.facets[i].arr[j].coordinates[k]
                        != actual.facets[i].arr[j].coordinates[k];
            }
        }
    }
    return differences;
}

void Check(const std::string &name, Particle &particle)
{
    particle.Rotate(0.31, 0.47, 0.19);
    particle.Resize(particle.MaximalDimention() * 1.07);

    const std::string path = "/tmp/mbs_raw_roundtrip_" + name + ".dat";
    particle.Output(path);

    Particle withSidecar;
    withSidecar.SetFromFile(path);
    const Differences sidecar = CountDifferences(particle, withSidecar);

    std::remove((path + ".visibility").c_str());
    Particle legacy;
    legacy.SetFromFile(path);
    const Differences legacyDifferences = CountDifferences(particle, legacy);

    std::remove(path.c_str());
    std::cout << name << ',' << particle.nFacets << ','
              << sidecar.topology << ',' << sidecar.coordinates << ',' << sidecar.visibility << ','
              << legacyDifferences.topology << ',' << legacyDifferences.coordinates << ','
              << legacyDifferences.visibility << std::endl;
}
}

int main()
{
    const complex ri(1.3116, 0.0);
    std::cout << "class,facets,sidecar_topology,sidecar_coordinates,sidecar_visibility,"
                 "legacy_topology,legacy_coordinates,legacy_visibility" << std::endl;
    { Hexagonal p(ri, 2.8, 4.0); Check("Hexagonal", p); }
    { ConcaveHexagonal p(ri, 2.8, 4.0, 35.0); Check("ConcaveHexagonal", p); }
    { Bullet p(ri, 2.8, 4.0, 1.0); Check("Bullet", p); }
    { BulletRosette p(ri, 1.0, 2.0, 0.5); Check("BulletRosette", p); }
    { Droxtal p(ri, DegToRad(32.35), DegToRad(71.81), 2.0); Check("Droxtal", p); }
    { HexagonalAggregate p(ri, 2.8, 4.0, 2); Check("HexagonalAggregate", p); }
    { CertainAggregate p(ri, 1.0); Check("CertainAggregate", p); }
    return 0;
}
