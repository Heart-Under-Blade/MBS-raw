#include "TracerGO.h"
#include "HandlerGO.h"
#include <iostream>

using namespace std;

TracerGO::TracerGO(Particle *particle, int reflNum, const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{
}

void TracerGO::TraceRandom(const AngleRange &betaRange, const AngleRange &gammaRange)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incomingEnergy = 0;
	m_outcomingEnergy = 0;
#endif

	vector<Beam> outBeams;
	double beta, gamma;

	CalcTimer timer;
	OutputStartTime(timer);

    long long orNum = gammaRange.number * betaRange.number;
    int betaNorm = (m_symmetry.beta < M_PI_2+FLT_EPSILON && m_symmetry.beta > M_PI_2-FLT_EPSILON) ? 1 : 2;
    double normIndex = gammaRange.number * betaNorm;
    // double norm = CalcNorm(orNum);
    m_handler->SetNormIndex(normIndex);

    double cs_beta = 0.0;
    long long count = 0;

    for (int i = 0; i <= betaRange.number; ++i)
	{
        beta = beta = betaRange.min + i*betaRange.step;
        CalcCsBeta(betaNorm, beta, betaRange, gammaRange, normIndex, cs_beta);

        for (int j = 0; j <= gammaRange.number; ++j)
		{
			gamma = (j + 0.5)*gammaRange.step;

            m_particle->Rotate(beta, gamma, 0);
			m_scattering->ScatterLight(beta, gamma, outBeams);
            m_handler->HandleBeams(outBeams, cs_beta);

#ifdef _CHECK_ENERGY_BALANCE
            m_incomingEnergy += m_scattering->GetIncedentEnergy()*cs_beta;
#endif
//			m_handler->WriteLog(to_string(i) + ", " + to_string(j) + " ");
//			OutputOrientationToLog(i, j, logfile);
            if (m_logTime == 0)
            {
                OutputProgress(orNum, ++count,
                               std::lround(RadToDeg(beta)),
                               std::lround(RadToDeg(gamma)), timer,
                               outBeams.size());
            }
            else
            {
                OutputProgress(orNum, ++count, i,j, timer, outBeams.size());
            }

            outBeams.clear();
		}

	}

    // m_incomingEnergy *= normIndex;
    m_handler->m_outputEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();
    m_handler->WriteMatricesToFile(m_resultDirName, 1000);
    OutputSummary(orNum, timer);
}

void TracerGO::TraceFixed(const double &beta, const double &gamma)
{
	double b = DegToRad(beta);
	double g = DegToRad(gamma);

	vector<Beam> outBeams;
	m_scattering->ScatterLight(b, g, outBeams);
    m_handler->HandleBeams(outBeams, 0);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

    m_handler->WriteMatricesToFile(m_resultDirName, 1000);
    //	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

void TracerGO::TraceMonteCarlo(const AngleRange &betaRange, const AngleRange &gammaRange, int nOrientations)
{
#ifdef _CHECK_ENERGY_BALANCE
    m_incomingEnergy = 0;
    m_outcomingEnergy = 0;
#endif
    vector<Beam> outBeams;
    double beta, gamma;

    CalcTimer timer;
    OutputStartTime(timer);

    string dir = CreateFolder(m_resultDirName);
    string fulldir = dir + m_resultDirName + '\\';

    long long nTacts;
    asm("rdtsc" : "=A"(nTacts));
    srand(nTacts);
    cout << endl << "NTacts = " << nTacts << endl;
//    srand(static_cast<unsigned>(time(0)));
    long long count = 0;

    for (int i = 0; i < nOrientations; ++i)
    {
        beta = RandomDouble(0, 1)*betaRange.max;
        gamma = RandomDouble(0, 1)*gammaRange.max;

        try
        {
            m_particle->Rotate(beta, gamma, 0);
            m_scattering->ScatterLight(beta, gamma, outBeams);
            m_handler->HandleBeams(outBeams, sin(beta));
#ifdef _CHECK_ENERGY_BALANCE
            m_incomingEnergy += m_scattering->GetIncedentEnergy()*sin(beta);
#endif
        }
        catch (...)
        {
        }

        outBeams.clear();

        OutputProgress(nOrientations, ++count, i, i, timer, outBeams.size());
    }

    double norm = CalcNorm(nOrientations);
    m_handler->SetNormIndex(norm);

    m_outcomingEnergy = ((HandlerGO*)m_handler)->ComputeTotalScatteringEnergy();

    m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;

    m_handler->WriteMatricesToFile(m_resultDirName, m_incomingEnergy);
    OutputSummary(nOrientations, timer);

}

double TracerGO::CalcNorm(long long orNum)
{
	double &symBeta = m_symmetry.beta;
    double tmp = /*(true) ? symBeta :*/ 1.0;
	double dBeta = -(cos(symBeta) - cos(0));
	return tmp/(orNum*dBeta);
}

void TracerGO::OutputSummary(int orNumber, CalcTimer &timer)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_summary += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber)
            /*+ "\nTotal scattering energy = " + to_string(D_tot)*/;

#ifdef _CHECK_ENERGY_BALANCE
    double passedEnergy = (m_handler->m_outputEnergy/m_incomingEnergy)*100;

    m_summary += "\nTotal incoming energy = " + to_string(m_incomingEnergy)
            + "\nTotal outcoming energy = " + to_string(m_handler->m_outputEnergy)
                 + " (S/4 = " + to_string(m_particle->Area()/4)
            + ")\nEnergy passed = " + to_string(passedEnergy) + '%';
#endif

	// out << "\nAveraged cross section = " << incomingEnergy*NRM;
    ofstream out(m_resultDirName+"_out.txt", ios::out);
	out << m_summary;
	out.close();

	cout << m_summary;
}
