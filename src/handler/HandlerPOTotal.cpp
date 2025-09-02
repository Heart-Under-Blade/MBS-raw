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
    outFile << "Theta 2pi*dcos M11 M12 M13 M14 "\
        "M21/ M22 M23 M24 "\
        "M31/ M32 M33 M34 "\
        "M41/ M42 M43 M44";
    matrix sum(4, 4);

    auto &Lp = *m_Lp;
    auto &Ln = *m_Ln;

    int &nT = m_sphere.nZenith;
    double &dT = m_sphere.zenithStep;

    int &nP = m_sphere.nAzimuth;

    for (int t = 0; t <= nT; ++t)
//    for (int t = nT; t >= 0; --t)
    {
        sum.Fill(0.0);
//        double tt = RadToDeg((m_sphere.zenithEnd - m_sphere.zenithStart) - t*dT);
        double tt = RadToDeg(m_sphere.zenithStart + (t*dT));

        for (int p = 0; p <= nP; ++p)
        {
            double radPhi = -p*m_sphere.azinuthStep;
            matrix m = M(p, t);
//#ifdef _DEBUG // DEB
//            if (t == nT && p == nP)
//                int fffff = 0;
//            double fff[4];
//            fff[0] = m[0][0];
//            fff[1] = m[0][1];
//            fff[2] = m[1][0];
//            fff[3] = m[1][1];
//#endif
            Lp[1][1] = cos(2*radPhi);
            Lp[1][2] = sin(2*radPhi);
            Lp[2][1] = -Lp[1][2];
            Lp[2][2] = Lp[1][1];

            Ln[1][2] = -Lp[1][2];
            Ln[2][1] = -Lp[2][1];

            if (t == 0)
            {
                sum += Lp*m*Lp;
            }
            else if (t == nT)
            {
                sum += Ln*m*Lp; // OPT: вынести Ln в отдельный случай
            }
            else
            {
                sum += m*Lp;
            }
        }

        double dS2 = (t == 0 || t == (nT)) ? 1.0-cos(0.5*dT)
                                           : (cos((t-0.5)*dT)-cos((t+0.5)*dT));

        sum /= m_sphere.nAzimuth;
#ifdef _DEBUG
        double dde = sum[0][0];
        std::cout << "@@@@@@@@@@@ " << dde << std::endl;
#endif
        dS2 *= M_2PI;
        outFile << std::endl << tt << ' ' << dS2 << ' ' /*<< nrg << ' '*/;
        outFile << sum;
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
