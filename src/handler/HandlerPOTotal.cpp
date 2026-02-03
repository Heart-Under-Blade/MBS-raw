#include "HandlerPOTotal.h"

#include "Mueller.hpp"
#include <iostream>
#include <iomanip>

HandlerPOTotal::HandlerPOTotal(Particle *particle, Light *incidentLight, int nTheta,
                               double wavelength)
    : HandlerPO(particle, incidentLight, nTheta, wavelength)
{
    betaMueller = new matrix(4, 4);
}

void HandlerPOTotal::WriteMatricesToFile(std::string &destName, double nrg)
{
    std::ofstream outFile(destName + ".dat", std::ios::out);

    if (!outFile.is_open())
    {
        // int fff = 0;
    }

    outFile << std::setprecision(10);
    outFile << "ScAngle 2pi*dcos "\
            "M11 M12 M13 M14 "\
            "M21 M22 M23 M24 "\
            "M31 M32 M33 M34 "\
            "M41 M42 M43 M44";
    matrix Msum(4, 4);

    auto &Lp = *m_Lp;
    auto &Ln = *m_Ln;

    int &nZen = m_sphere.nZenith;
    double &dZen = m_sphere.zenithStep;
    int &nAz = m_sphere.nAzimuth;

    for (int iZen = 0; iZen <= nZen; ++iZen)
//    for (int iZen = nZen; iZen >= 0; --iZen)
    {
        Msum.Fill(0.0);
        double radZen = m_sphere.zenithStart + iZen*dZen;
//        double tt = RadToDeg(m_sphere.zenithEnd) - RadToDeg((t*dT));

        for (int iAz = 0; iAz <= nAz; ++iAz)
        {
            double radAz = -iAz*m_sphere.azinuthStep;
            matrix m = M(iAz, iZen);

            Lp[1][1] = cos(2*radAz);
            Lp[1][2] = sin(2*radAz);
            Lp[2][1] = -Lp[1][2];
            Lp[2][2] = Lp[1][1];

            Ln = Lp;
            Ln[1][2] *= -1;
            Ln[2][1] *= -1;

            if (radZen > M_PI-__FLT_EPSILON__)
            {
                Msum += Lp*m*Lp;
            }
            else if (radZen < __FLT_EPSILON__)
            {
                Msum += Ln*m*Lp;
            }
            else
            {
                Msum += m*Lp;
            }
        }

        double _2Pi_dcos = (radZen > M_PI-__FLT_EPSILON__ || radZen < __FLT_EPSILON__) ?
                    1.0-cos(0.5*dZen) :
                    cos((iZen-0.5)*dZen)-cos((iZen+0.5)*dZen);

        Msum /= m_sphere.nAzimuth;
        _2Pi_dcos *= M_2PI;
        outFile << std::endl << RadToDeg(radZen) << ' ' << _2Pi_dcos << ' ' /*<< nrg << ' '*/;
        outFile << Msum;
    }

    outFile.close();
}

void HandlerPOTotal::AddToMueller()
{
#ifdef _DEBUG // DEB
    double sum = 0;
#endif
    for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
    {
        auto &diffM = m_diffractedMatrices[q];

        for (int t = 0; t <= m_sphere.nZenith; ++t)
        {
            for (int p = 0; p <= m_sphere.nAzimuth; ++p)
            {
                matrix m = Mueller(diffM(p, t));
#ifdef _DEBUG // DEB
                double fff[4];
                fff[0] = m[0][0];
#endif
                if (q == 0 && t == 0 && p == 0)
                {
                    *betaMueller += m*normIndexGamma;
                }

                m *= m_sinZenith;

#ifdef _DEBUG // DEB
                complex ddd[4];
                ddd[0] = diffM(p, t)[0][0];
                double d = m[0][0];
//                if (t == 160)
                {
                    sum += d;
//                    m_logFile << p << ' ' << t << ' ' << sum << std::endl;
                }
#endif
                M.insert(p, t, m);
            }
        }
    }
#ifdef _DEBUG // DEB
    double fffefwe = M(0,0)[0][0];
    int ff = 0;
#endif
}

void HandlerPOTotal::OutputContribution(double angle, double energy)
{
    *(betaFile) << RadToDeg(angle) << ' ' << energy << ' ';
    *(betaFile) << *(betaMueller) << std::endl;
    betaMueller->Fill(0.f);
}
