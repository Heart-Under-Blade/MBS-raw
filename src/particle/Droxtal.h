#pragma once

#include "Hexagonal.h"

class Droxtal : public Hexagonal
{
public:
	Droxtal(const complex &refrIndex,
			double topAngle, double midAngle, double radius);
};
