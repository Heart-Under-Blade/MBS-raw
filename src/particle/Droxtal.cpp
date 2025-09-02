#include "Droxtal.h"

Droxtal::Droxtal(const complex &refrIndex,
				 double topAngle, double midAngle, double radius)
{
	isConcave = false;

	const double halfHeight = radius*cos(topAngle);
	const double midHalfHeight = radius*cos(midAngle);

	const double topRadius = radius*sin(topAngle);
	const double midRadius = radius*sin(midAngle);

	SetSize(midRadius*2, halfHeight*2);
	Init(20, refrIndex);

	defaultFacets[0] = CreateBase(6, topRadius, halfHeight);
	defaultFacets[19] = CreateBase(6, topRadius, -halfHeight);

	Facet topMidBase = CreateBase(6, midRadius, midHalfHeight);
	Facet botMidBase = CreateBase(6, midRadius, -midHalfHeight);

	int iFacet = 0;

	CreateSideFacets(defaultFacets[0], topMidBase, iFacet);
	CreateSideFacets(topMidBase, botMidBase, iFacet);
	CreateSideFacets(botMidBase, defaultFacets[19], iFacet);

	SetDefaultNormals();
	Reset();
	SetDefaultCenters();
}
