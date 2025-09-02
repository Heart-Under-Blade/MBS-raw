#include "Facet.h"
#include "Intersection.h"

void Facet::SetNormal()
{
	ex_normal = Normal();
	in_normal = -ex_normal;
}

void Facet::SetCenter()
{
    center = Center();
}

bool Facet::IsConormal(Point3f normal) const
{
    __m128 dp = _dp(_set(this->normal[Out]), _set(normal), MASK_FULL);
    return dp[0] > __FLT_EPSILON__;
}

Facet &Facet::operator =(const Facet &other)
{
	if (this != &other)
	{
		Polygon::operator =(other);
		in_normal = other.in_normal;
		ex_normal = other.ex_normal;
		center = other.center;
	}

	return *this;
}

