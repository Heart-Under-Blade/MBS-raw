#include "Particle.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include "global.h"

using namespace std;

Particle::Particle()
{
    isConcave = false;
}

void Particle::SetFromFile(const std::string &filename)
{
    std::ifstream pfile(filename, std::ios::in);

    if (!pfile.is_open())
    {
        std::cerr << "File \"" << filename << "\" is not found" << std::endl;
        throw std::exception();
    }

    const int bufSize = 1024;
    char *buff = (char*)malloc(sizeof(char) * bufSize);

    nFacets = 0;
    Facet *facet = &(defaultFacets[nFacets++]);

    char *ptr, *trash;

    pfile.getline(buff, bufSize);
    isConcave = strtol(buff, &trash, 10);

    pfile.getline(buff, bufSize);
    ptr = strtok(buff, " ");
    isAggregated = strtol(ptr, &trash, 10);

    if (isAggregated)
    {
        ptr = strtok(NULL, " ");

        if (ptr != nullptr)
        {
            nFacetsInPart = strtod(ptr, &trash);
        }
    }

    // read symmetry params
    {
        pfile.getline(buff, bufSize);

        ptr = strtok(buff, " ");
        m_symmetry.beta = DegToRad(strtod(ptr, &trash));

        ptr = strtok(NULL, " ");
        m_symmetry.gamma = DegToRad(strtod(ptr, &trash));
    }

    pfile.getline(buff, bufSize); // skip empty line

    while (!pfile.eof()) // read vertices of facets
    {
        pfile.getline(buff, bufSize);
        ptr = strtok(buff, " ");

        if (strlen(buff) == 0 || !strcmp(ptr, "\r"))
        {
            facet = &(defaultFacets[nFacets++]);
            continue;
        }

        int c_i = 0;

        while (ptr != NULL)
        {
            facet->arr[facet->nVertices].coordinates[c_i++] = strtod(ptr, &trash);
            ptr = strtok(NULL, " ");
        }

        ++(facet->nVertices);
    }

    pfile.close();

    // correction of number of facet
    while (defaultFacets[nFacets-1].nVertices == 0)
    {
        --nFacets;
    }

    if (isConcave)
    {
//        RemoveWalls();
    }

    SetDefaultNormals();
    Reset();
    SetDefaultCenters();

    if (isConcave || isAggregated)
    {
        for (int i = 0; i < nFacets; ++i)
        {
            defaultFacets[i].isVisibleIn = false;
            defaultFacets[i].isVisibleOut = false;
            facets[i].isVisibleIn = false;
            facets[i].isVisibleOut = false;
        }
    }
}

void Particle::Init(int facetCount, const complex &refrIndex)
{
    nFacets = facetCount;
    m_refractiveIndex = refrIndex;
}

void Particle::RotateCenters()
{
    for (int i = 0; i < nFacets; ++i)
    {
        RotatePoint(defaultFacets[i].center, facets[i].center);
    }
}

void Particle::Rotate(double beta, double gamma, double alpha)
{
    rotAngle = Angle{alpha, beta, gamma};
    SetRotateMatrix(beta, gamma, alpha);

    // REF: слить всё в один цикл
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            RotatePoint(defaultFacets[i].arr[j], facets[i].arr[j]);
        }
    }

    RotateNormals();
    RotateCenters();
}

void Particle::Fix()
{
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            defaultFacets[i].arr[j] = facets[i].arr[j];
        }
    }
}

void Particle::Scale(double ratio)
{
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            defaultFacets[i].arr[j].cx *= ratio;
            defaultFacets[i].arr[j].cy *= ratio;
            defaultFacets[i].arr[j].cz *= ratio;
        }
    }
}

// vector<unsigned> Particle::FindFlippedFacets() const
// {
//     vector<unsigned> facets;

//     if (isAggregated)
//     {
//         vector<Particle> parts;
//         DisconnectParts(parts);

//         unsigned totalNFacets = 0;

//         for (unsigned i = 0; i < parts.size(); ++i)
//         {
//             vector<unsigned> partFacets;
//             parts[i].FindFlippedFacets(partFacets);

//             for (unsigned j = 0; j < partFacets.size(); ++j)
//             {
//                 partFacets[j] += totalNFacets;
//                 facets.push_back(partFacets[j]);
//             }

//             totalNFacets += parts[i].nFacets;
//         }
//     }
//     else
//     {
//         FindFlippedFacets(facets);
//     }

//     return facets;
// }

// void Particle::DisconnectParts(std::vector<Particle> &parts) const
// {
//     for (unsigned i = 0; i < partIndices.size(); ++i)
//     {
//         Particle part;

//         for (unsigned j = partIndices[i].from; j < partIndices[i].to; ++j)
//         {
//             part.AddFacet(facets[j]);
//         }

//         part.UpdateFacets();
//         parts.push_back(part);
//     }
// }

void Particle::Resize(double size)
{
    double ratio = (isAggregated) ? size/MaximalDimentionPart() :
                                    size/MaximalDimention();
    Scale(ratio);
    Reset();
}

Facet Particle::CreateBase(int nVertices, double radius, double zCoordinate)
{
    Facet base;
    double step = M_2PI/nVertices;

    int begin;
    int end;
    int inc;

    if (zCoordinate >= 0)
    {
        begin = 0;
        end = nVertices;
        inc = 1;
        base.normal[0] = Point3f(0, 0, -1);
        base.normal[1] = Point3f(0, 0, 1);
    }
    else  // reverse vertex order
    {
        begin = nVertices-1;
        end = -1;
        inc = -1;
        base.normal[0] = Point3f(0, 0, 1);
        base.normal[1] = Point3f(0, 0, -1);
    }

    for (int i = begin; i != end; i += inc)
    {
        double a = step*i;
        base.AddVertex(Point3f(radius*cos(a), radius*sin(a), zCoordinate));
    }

    return base;
}

void Particle::CreateSideFacets(const Facet &top, const Facet &bottom,
                                int &iFacet)
{
    std::vector<Facet> sideFacets;
    sideFacets = ConnectBases(top, bottom);

    for (size_t i = 0; i < sideFacets.size(); ++i)
    {
        defaultFacets[++iFacet] = sideFacets[i];
    }
}

std::vector<Facet> Particle::ConnectBases(const Facet &top,
                                          const Facet &bottom)
{
    std::vector<Facet> faces;

    if (top.nVertices == bottom.nVertices)
    {
        bool isInverse = DotProduct(bottom.normal[1], top.normal[1]) < 0;
        int end = top.nVertices-1;

        int i1 = end;
        int i2 = 0;

        if (top.normal[1].coordinates[2] >= 0)
        {
            for (int i = 0; i < top.nVertices; ++i)
            {
                Facet f;

                f.AddVertex(top.arr[i2]);
                f.AddVertex(top.arr[i1]);

                if (isInverse)
                {
                    f.AddVertex(bottom.arr[end-i1]);
                    f.AddVertex(bottom.arr[end-i2]);
                }
                else
                {
                    f.AddVertex(bottom.arr[i1]);
                    f.AddVertex(bottom.arr[i2]);
                }

                faces.push_back(f);

                i1 = i2;
                ++i2;
            }
        }
        else
        {
            for (int i = 0; i < top.nVertices; ++i)
            {
                Facet f;

                f.AddVertex(top.arr[i1]);
                f.AddVertex(top.arr[i2]);

                if (isInverse)
                {
                    f.AddVertex(bottom.arr[end-i2]);
                    f.AddVertex(bottom.arr[end-i1]);
                }
                else
                {
                    f.AddVertex(bottom.arr[i2]);
                    f.AddVertex(bottom.arr[i1]);
                }

                faces.push_back(f);

                i1 = i2;
                ++i2;
            }
        }
    }
    else
    {
        throw "Number of base vertices are not equal";
    }

    return faces;
}

void Particle::FindFlippedFacets(std::vector<unsigned> &facetIndices) const
{
    Point3f globalCenter = Center();

    for (unsigned i = 0; i < nFacets; ++i)
    {
        if (!facets[i].IsConormal(facets[i].center - globalCenter))
        {
            facetIndices.push_back(i);
        }
    }
}

void Particle::Concate(const std::vector<Particle> &parts)
{
    int i = 0;
    nFacets = 0;

    for (const Particle &part : parts)
    {
        nFacets += part.nFacets;

        for (int j = 0; j < part.nFacets; ++j)
        {
            defaultFacets[i++] = part.facets[j];
        }
    }

    isAggregated = true;
}

double Particle::LongRadius() const
{
    Point3f p0(0, 0, 0);

    double radius = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            Point3f v_len = facets[i].arr[j] - p0;
            double len = Length(v_len);

            if (len > radius)
            {
                radius = len;
            }
        }
    }

    return radius;
}

double Particle::MaximalDimention() const
{
    double Dmax = 0;
    double newDmax;

    std::vector<Point3f> vertices;

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            vertices.push_back(defaultFacets[i].arr[j]);
        }
    }

    for (int i = 0; i < vertices.size(); ++i)
    {
        for (int j = 0; j < vertices.size(); ++j)
        {
            newDmax = Length(vertices[j] - vertices[i]);

            if (newDmax > Dmax)
            {
                Dmax = newDmax;
            }
        }
    }

    return Dmax;
}

double Particle::MaximalDimentionPart() const
{
    double Dmax = 0;
    double newDmax;

    Polygon512 pol;

    for (int i = 0; i < nFacetsInPart; ++i)
    {
        pol.Concat(defaultFacets[i]);
    }

    for (int i = 0; i < pol.nVertices; ++i)
    {
        for (int j = 0; j < pol.nVertices; ++j)
        {
            newDmax = Length(pol.arr[j] - pol.arr[i]);

            if (newDmax > Dmax)
            {
                Dmax = newDmax;
            }
        }
    }

    return Dmax;
}

double Particle::Area()
{
    double area = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        area += facets[i].Area();
    }

    return area;
}

Point3f Particle::Center() const
{
    Point3f center = Point3f(0, 0, 0);
    int nVertices = 0;

    for (int i = 0; i < nFacets; ++i)
    {
        auto &facet = defaultFacets[i];
        nVertices += facet.nVertices;

        for (int j = 0; j < facet.nVertices; ++j)
        {
            center = center + facet.arr[j];
        }
    }

    center = center/nVertices;
    return center;
}

double Particle::Volume()
{
    double volume = 0;
    Point3f center = Center();

    for (int i = 0; i < nFacets; ++i)
    {
        const Facet &facet = defaultFacets[i];

        Point3f p = ProjectPointToPlane(center, facet.ex_normal,
                                                  facet.in_normal);
        double h = Length(p - center);
        volume += (facet.Area()*h)/3;
    }

    return volume;
}

const complex &Particle::GetRefractiveIndex() const
{
    return m_refractiveIndex;
}

const Symmetry &Particle::GetSymmetry() const
{
    return m_symmetry;
}

void Particle::Move(float dx, float dy, float dz)
{
    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < defaultFacets[i].nVertices; ++j)
        {
            facets[i].arr[j] = defaultFacets[i].arr[j] + Point3f(dx, dy, dz);
        }
    }
}

bool Particle::IsConcave() const
{
    return isConcave;
}

void Particle::Output(std::string name)
{
    std::ofstream file(name, std::ios::out);

    file << (int)isConcave << std::endl;
    file << (int)isAggregated;

    if (isAggregated)
    {
        file << ' ' << nFacetsInPart;
    }

    file << std::endl;

    file << RadToDeg(m_symmetry.beta) << ' '
           << RadToDeg(m_symmetry.gamma)
           << std::endl << std::endl;

    for (int i = 0; i < nFacets; ++i)
    {
        for (int j = 0; j < facets[i].nVertices; ++j)
        {
            Point3f p = facets[i].arr[j];
            file << p.coordinates[0] << ' '
                            << p.coordinates[1] << ' '
                            << p.coordinates[2] << ' '
                            /*<< i */;
            file << std::endl;
        }

        if ((i + 1) != nFacets)
        {
            file << std::endl;
        }
    }

    file.close();
}

void Particle::SetRefractiveIndex(const complex &value)
{
    m_refractiveIndex = value;
}

void Particle::SetDefaultNormals()
{
    for (int i = 0; i < nFacets; ++i)
    {
        defaultFacets[i].SetNormal();
    }
}

void Particle::Reset()
{
    for (int i = 0; i < nFacets; ++i)
    {
        facets[i] = defaultFacets[i];
    }
}

void Particle::SetDefaultCenters()
{
    for (int i = 0; i < nFacets; ++i)
    {
        defaultFacets[i].SetCenter();
    }
}

void Particle::SetRotateMatrix(double beta, double gamma, double alpha)
{
    double cosA, cosB, cosG,
            sinA, sinB, sinG;

    sincos(alpha, &sinA, &cosA);
    sincos(beta,  &sinB, &cosB);
    sincos(gamma, &sinG, &cosG);

    double cosAcosB = cosA*cosB;
    double sinAcosG = sinA*cosG;
    double sinAsinG = sinA*sinG;

    m_rotMatrix[0][0] = cosAcosB*cosG - sinAsinG;
    m_rotMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
    m_rotMatrix[2][0] = -sinB*cosG;

    m_rotMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
    m_rotMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
    m_rotMatrix[2][1] = sinB*sinG;

    m_rotMatrix[0][2] = cosA*sinB;
    m_rotMatrix[1][2] = sinA*sinB;
    m_rotMatrix[2][2] = cosB;
}

void Particle::RemoveWalls()
{
    for (int i = 0; i < nFacets; ++i)
    {
        auto &f1 = defaultFacets[i];

        for (int j = i+1; j < nFacets; ++j)
        {
            auto &f2 = defaultFacets[j];

            bool isFoundF = true;

            if (f1.nVertices == f2.nVertices)
            {
                for (int m = 0; m < f1.nVertices; ++m)
                {
                    bool isFoundV = false;

                    for (int k = 0; k < f2.nVertices; ++k)
                    {
                        if (f1.arr[m].IsEqualTo(f2.arr[k], 0.005))
                        {
                            isFoundV = true;
                            break;
                        }
                    }

                    if (!isFoundV)
                    {
                        isFoundF = false;
                        break;
                    }
                }

                if (isFoundF)
                {
                    defaultFacets[i].nVertices = 0;
                    defaultFacets[j].nVertices = 0;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < nFacets; ++i)
    {
        if (defaultFacets[i].nVertices == 0)
        {
            for (int j = i+1; j < nFacets; ++j)
            {
                defaultFacets[j-1] = defaultFacets[j];
            }

            --nFacets;
            --i;
        }
    }
}

void Particle::RotateNormals()
{
    for (int i = 0; i < nFacets; ++i)
    {
        RotatePoint(defaultFacets[i].in_normal, facets[i].in_normal);
    }

    SetDParams();

    for (int i = 0; i < nFacets; ++i)
    {
        facets[i].ex_normal = -facets[i].in_normal;
        facets[i].ex_normal.d_param = -facets[i].in_normal.d_param;
    }
}

void Particle::SetDParams()
{
    for (int i = 0; i < nFacets; ++i)
    {
        double d = DotProduct(facets[i].arr[0], facets[i].in_normal);
        facets[i].in_normal.d_param = -d;
    }
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
    result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
    result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
    result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}

void Particle::SetSymmetry(double beta, double gamma, double alpha)
{
    m_symmetry.beta = beta;
    m_symmetry.gamma = gamma;
    m_symmetry.alpha = alpha;
}
