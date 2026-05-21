#include "TracerPOTotal.h"
#include "HandlerPOTotal.h"

#include <iostream>

using namespace std;

TracerPOTotal::TracerPOTotal(Particle *particle, int nActs,
                             const string &resultFileName)
    : TracerPO(particle, nActs, resultFileName)
{
}

void TracerPOTotal::TraceRandom(const AngleRange &betaRange,
                                const AngleRange &gammaRange)
{
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy = 0;
    m_outcomingEnergy = 0;
#endif

    CalcTimer timer;
    long long count = 0;
    long long nOrientations = (betaRange.number+1) * (gammaRange.number);

    vector<Beam> outBeams;
    double beta, gamma;

    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normIndex = gammaRange.number * betaNorm;
    m_handler->SetNormIndex(normIndex);

    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + "/" + m_resultDirName;
#endif

    m_handler->betaFile = new std::ofstream(m_resultDirName + "_beta.dat", ios::out);

    *(m_handler->betaFile) << "Beta CS M11 M12 M13 M14 "\
        "M21 M22 M23 M24 "\
        "M31 M32 M33 M34 "\
        "M41 M42 M43 M44";

    timer.Start();
    OutputStartTime(timer);

    for (int i = 0; i <= betaRange.number; ++i)
    {
        beta = betaRange.min + i*betaRange.step;

        double dcos;
        bool isPole = CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normIndex, dcos);
        m_handler->SetSinZenith(dcos);

        for (int j = 0; j < gammaRange.number; ++j)
        {
//            if (j != 0 && isPole)
//            {
//                continue;
//            }

            gamma = gammaRange.min + (j+0.5)*gammaRange.step;
            m_particle->Rotate(/*M_PI-*/beta, /*M_PI+*/gamma, 0);

            if (!shadowOff)
            {
                m_scattering->FormShadowBeam(outBeams);
            }

            bool ok = m_scattering->ScatterLight(0, 0, outBeams);

            if (ok)
            {
                m_handler->HandleBeams(outBeams, dcos);
            }
            else
            {
                std::cout << std::endl << "Orientation (" << i << ", " << j << ") has been skipped!!!" << std::endl;
            }

            m_incomingEnergy += m_scattering->GetIncedentEnergy()*dcos;

            OutputProgress(dir, nOrientations, count, i, j, timer, outBeams.size());
            outBeams.clear();
            ++count;
        }

        static_cast<HandlerPOTotal*>(m_handler)->OutputContribution(beta, m_incomingEnergy);
    }

    static_cast<HandlerPOTotal*>(m_handler)->betaFile->close();


    EraseConsoleLine(60);
    std::cout << "100%" << std::endl;

    // m_handler->m_outputEnergy = m_handler->ComputeTotalScatteringEnergy();
    m_handler->WriteTotalMatricesToFile(m_resultDirName);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy, true);
    string incohName = m_resultDirName+string("_incoh");
    m_handler->WriteMatricesToFile(incohName, m_incomingEnergy, false);

    OutputStatisticsPO(timer, nOrientations, m_resultDirName);
}

void TracerPOTotal::TraceMonteCarlo(const AngleRange &betaRange,
                                    const AngleRange &gammaRange,
                                    int nOrientations)
{
    CalcTimer timer;
    long long count = 0;

    ofstream outFile(m_resultDirName + ".dat", ios::out);

    if (!outFile.is_open())
    {
        std::cerr << "Error! File \"" << m_resultDirName << "\" was not opened. "
                  << __FUNCTION__;

        throw std::exception();
    }

    vector<Beam> outBeams;
    double beta, gamma;
    timer.Start();

    long long nTacts;
    asm("rdtsc" : "=A"(nTacts));
    srand(nTacts);
//    srand(static_cast<unsigned>(time(0)));

    for (int i = 0; i < nOrientations; ++i)
    {
        beta = RandomDouble(0, 1)*betaRange.max;
        gamma = RandomDouble(0, 1)*gammaRange.max;

        m_particle->Rotate(beta, gamma, 0);
        m_scattering->ScatterLight(beta, gamma, outBeams);

        m_handler->HandleBeams(outBeams, sin(beta));
        outBeams.clear();

        ++count;
        OutputProgress(m_resultDirName, nOrientations, count, i, i, timer, outBeams.size());
    }

    m_handler->WriteTotalMatricesToFile(m_resultDirName);

    std::string dir = CreateFolder(m_resultDirName);
    m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy, true);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy, false);
    OutputStatisticsPO(timer, nOrientations, dir);
    outFile.close();
}
