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
    long long nOrientations = (betaRange.number) * (gammaRange.number);

#ifdef _DEBUG  /* DEB */
    ofstream outFile(m_resultDirName + "log.dat", ios::out);
    if (!outFile.is_open())
    {
        std::cerr << "Error! File \"" << m_resultDirName
                  << "\" was not opened. " << __FUNCTION__;
        throw std::exception();
    }
#endif


    vector<Beam> outBeams;
    double beta, gamma;

    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normIndex = gammaRange.number * betaNorm;
    m_handler->SetNormIndex(normIndex);

    std::string dir = CreateFolder(m_resultDirName);
#ifdef _WIN32
    m_resultDirName += '\\' + m_resultDirName;
#else
    m_resultDirName = dir + m_resultDirName;
#endif

    m_handler->betaFile = new std::ofstream(m_resultDirName + "_beta.dat", ios::out);

    *(m_handler->betaFile) << "Beta CS M11 M12 M13 M14 "\
        "M21 M22 M23 M24 "\
        "M31 M32 M33 M34 "\
        "M41 M42 M43 M44";

//    ++nOrientations;
    timer.Start();
    OutputStartTime(timer);

#ifdef _DEBUG  /* DEB */
    for (int i = 0; i <= betaRange.number; ++i)
    {
//        std::cout  << "i: "<< i << std::endl;
#else
    for (int i = 0; i <= betaRange.number; ++i)
    {
#endif
        beta = betaRange.min + i*betaRange.step;

//		double dcos = (i == 0 || i == betaRange.number)
//				? (1.0-cos(0.5*betaRange.step))/normIndex
//				: (cos((i-0.5)*betaRange.step) -
//				   cos((i+0.5)*betaRange.step))/normIndex;

        double dcos;
        CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normIndex, dcos);
        m_handler->SetSinZenith(dcos);

#ifdef _DEBUG // DEB
        for (int j = 0; j < gammaRange.number; ++j)
        {
//            std::cout  << "j: "<< j << std::endl;
//			beta = DegToRad(15); gamma = DegToRad(0);
#else
        for (int j = 0; j < gammaRange.number; ++j)
        {
#endif
            gamma = gammaRange.min + j*gammaRange.step;
            m_particle->Rotate(/*M_PI-*/beta, /*M_PI+*/gamma, 0);

            if (!shadowOff)
            {
                m_scattering->FormShadowBeam(outBeams);
            }

#ifdef _DEBUG  /* DEB */
            // if (i == 4 && j ==0)
            //     int ffff = 0;
#endif
            bool ok = m_scattering->ScatterLight(0, 0, outBeams);

            if (ok)
            {
#ifdef _DEBUG  /* DEB */
//             vector<Beam> be = outBeams;
//            for (int k = 1; k < be.size(); ++k) {
//                if (be[k].nActs==2) {
//                    outBeams.push_back(be[k]);
//                }
//            }
#endif
                m_handler->HandleBeams(outBeams, dcos);
//#ifdef _DEBUG  /* DEB */
//            double sum = 0;
//            for (int k = 0; k < outBeams.size(); ++k) {
//                double aaa = outBeams[k].Area();
////                ofstream outFile(m_resultDirName + to_string(k) + "_vertices.dat", ios::out);
//                sum += aaa;
//                std::vector<int> tr;
//                m_handler->m_tracks->RecoverTrack(outBeams[k], m_particle->nFacets, tr);
//                outFile << outBeams[k].id << " ";
//                for (int l = 0; l < tr.size(); ++l) {
//                    outFile << tr[l] << " ";
//                }
//               outFile << aaa << std::endl;
////                outFile << outBeams[k].arr[0] << endl << endl;
////                outFile.close();
//            }
//            outFile.close();
//            double mm = ((HandlerPO*)m_handler)->M(0, 0)[0][0];
//            std::cout << mm << " " << j << std::endl;
//#else
//#endif
            }
            else
            {
                std::cout << std::endl << "Orientation (" << i << ", " << j << ") has been skipped!!!" << std::endl;
            }

            m_incomingEnergy += m_scattering->GetIncedentEnergy()*dcos;

            OutputProgress(nOrientations, count, i, j, timer, outBeams.size());
            outBeams.clear();
//            std::cout << "sdfewtwr";
            ++count;
        }

        static_cast<HandlerPOTotal*>(m_handler)->OutputContribution(beta, m_incomingEnergy);
    }

    static_cast<HandlerPOTotal*>(m_handler)->betaFile->close();

//#ifdef _DEBUG  /* DEB */
//    double mm = ((HandlerPO*)m_handler)->M(0, 0)[0][0];
//    outFile.close();
//#endif

    EraseConsoleLine(60);
    std::cout << "100%" << std::endl;

    // m_handler->m_outputEnergy = m_handler->ComputeTotalScatteringEnergy();
    m_handler->WriteTotalMatricesToFile(m_resultDirName);
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
//#ifndef _DEBUG

    OutputStatisticsPO(timer, nOrientations, m_resultDirName);
//#endif
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
        OutputProgress(nOrientations, count, i, i, timer, outBeams.size());
    }

    m_handler->WriteTotalMatricesToFile(m_resultDirName);

    std::string dir = CreateFolder(m_resultDirName);
    m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputStatisticsPO(timer, nOrientations, dir);
    outFile.close();
}
