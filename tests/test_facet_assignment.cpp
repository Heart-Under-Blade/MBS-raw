#include "geometry/Facet.h"
#include "particle/Bullet.h"
#include "particle/ConcaveHexagonal.h"
#include "particle/HexagonalAggregate.h"

#include <cassert>

namespace
{
void CheckConcaveHexagonal(const ConcaveHexagonal &particle)
{
    assert(particle.nFacets == 18);
    for (int i = 0; i < particle.nFacets; ++i)
    {
        const bool cavityFacet = i < 6 || i >= 12;
        const bool sideFacet = i >= 6 && i < 12;
        assert(particle.facets[i].isVisibleOut == !cavityFacet);
        assert(particle.facets[i].isVisibleIn == !sideFacet);
    }
}

void CheckBullet(const Bullet &particle)
{
    assert(particle.nFacets == 13);
    for (int i = 0; i < particle.nFacets; ++i)
    {
        assert(!particle.facets[i].isVisibleIn);
        assert(!particle.facets[i].isVisibleOut);
    }
}

void CheckHexagonalAggregate(const HexagonalAggregate &particle)
{
    const int hiddenOutside[] = {1, 2, 3, 7, 10, 11, 12, 15};
    bool expectedOutside[16] = {};
    for (int i = 0; i < 16; ++i)
    {
        expectedOutside[i] = true;
    }
    for (unsigned i = 0; i < sizeof(hiddenOutside)/sizeof(hiddenOutside[0]); ++i)
    {
        expectedOutside[hiddenOutside[i]] = false;
    }

    assert(particle.nFacets == 16);
    for (int i = 0; i < particle.nFacets; ++i)
    {
        assert(!particle.facets[i].isVisibleIn);
        assert(particle.facets[i].isVisibleOut == expectedOutside[i]);
    }
}
}

int main()
{
    Facet source;
    source.isVisibleIn = false;
    source.isVisibleOut = false;

    Facet destination;
    assert(destination.isVisibleIn);
    assert(destination.isVisibleOut);

    destination = source;
    assert(destination.isVisibleIn == source.isVisibleIn);
    assert(destination.isVisibleOut == source.isVisibleOut);

    const complex refractiveIndex(1.31, 0.0);

    ConcaveHexagonal concave(refractiveIndex, 10.0, 20.0, 30.0);
    CheckConcaveHexagonal(concave);
    concave.Rotate(0.31, 0.47, 0.19);
    CheckConcaveHexagonal(concave);
    concave.Resize(1.2);
    CheckConcaveHexagonal(concave);

    Bullet bullet(refractiveIndex, 10.0, 20.0, 5.0);
    CheckBullet(bullet);
    bullet.Rotate(0.31, 0.47, 0.19);
    CheckBullet(bullet);
    bullet.Resize(1.2);
    CheckBullet(bullet);

    HexagonalAggregate aggregate(refractiveIndex, 10.0, 20.0, 2);
    CheckHexagonalAggregate(aggregate);
    aggregate.Rotate(0.31, 0.47, 0.19);
    CheckHexagonalAggregate(aggregate);
    aggregate.Resize(1.2);
    CheckHexagonalAggregate(aggregate);
    return 0;
}
