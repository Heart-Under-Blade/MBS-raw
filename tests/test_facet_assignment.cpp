#include "geometry/Facet.h"

#include <cassert>

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
    return 0;
}
