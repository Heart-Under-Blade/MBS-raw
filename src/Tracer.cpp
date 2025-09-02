#include "Tracer.h"

#include <ostream>
#include <string>
#include <iostream>
#include <assert.h>
#include "global.h"
#include "macro.h"
//#ifdef _TRACK_ALLOW
//ofstream trackMapFile("tracks_deb.dat", ios::out);
//#endif

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"

using namespace std;

Tracer::Tracer(Particle *particle, int nActs, const string &resultFileName)
    : m_resultDirName(resultFileName)
{
    SetIncidentLight(particle);

    if (particle->IsConcave())
    {
        m_scattering = new ScatteringNonConvex(particle, &m_incidentLight, true, nActs);
    }
    else
    {
        m_scattering = new ScatteringConvex(particle, &m_incidentLight, true, nActs);
    }

    m_particle = m_scattering->m_particle;
    m_symmetry = m_particle->GetSymmetry();
}

Tracer::~Tracer()
{
}

void Tracer::CalcCsBeta(int betaNorm, double beta, const AngleRange &betaRange,
                          const AngleRange &gammaRange, double normIndex, double &cs_beta)
{
    if (betaNorm == 2)
        cs_beta = (fabs(beta)<=FLT_EPSILON || fabs(beta-M_PI)<=FLT_EPSILON) ?
                      (1.0-cos(0.5*betaRange.step)) :
                      (cos(beta-0.5*betaRange.step)-cos(beta+0.5*betaRange.step));
    else
    {
        cs_beta = (fabs(beta)<=FLT_EPSILON) ?
                      (1.0-cos(0.5*betaRange.step)) :
                      (cos(beta-0.5*betaRange.step)-cos(beta+0.5*betaRange.step));
        if(fabs(beta-M_PI/2.0)<=FLT_EPSILON)
            cs_beta *= 0.5;
    }

    cs_beta = cs_beta/normIndex;

    if(fabs(beta)<=FLT_EPSILON || fabs(beta-M_PI)<=FLT_EPSILON)
        cs_beta *= gammaRange.number;
}

void Tracer::SetIncidentLight(Particle *particle)
{
    m_incidentLight.direction = Point3f(0, 0, -1);
    m_incidentLight.polarizationBasis = Point3f(0, 1, 0);

    Point3f point = m_incidentLight.direction * particle->LongRadius();
    m_incidentLight.direction.d_param = DotProduct(point, m_incidentLight.direction);
}

void Tracer::OutputOrientationToLog(int i, int j, ostream &logfile)
{
    logfile << "i: " << i << ", j: " << j << endl;
    logfile.flush();
}

void Tracer::OutputProgress(int nOrientation, long long count,
                            int zenith, int azimuth, CalcTimer &timer, int nBeams)
{
    string split = "\t";

    auto now = timer.SecondsElapsed();

    if (now - m_timeElapsed > m_logTime)
    {
        m_timeElapsed = now;
//        EraseConsoleLine(60);

        string progressLine = to_string((count*100)/nOrientation) +
                "%" + split + "(" + to_string(zenith) + ";" +
                to_string(azimuth) + ")" + split + timer.Elapsed() + split +
                m_resultDirName+split+"beams:"+split+to_string(nBeams);

        cout << progressLine;
        RenameConsole(progressLine);
        cout << endl;
    }
}


double Tracer::CalcNorm(long long orNum)
{
//	return 1.0/(double)orNum;

    double &symBeta = m_symmetry.beta;

    if (symBeta < M_PI - FLT_EPSILON) //
    {
        double tmp = symBeta;
        double dBeta = -(cos(symBeta) - 1);
        dBeta = (dBeta < 0) ? -dBeta : dBeta;
        return tmp/(orNum*dBeta);
    }
    else // otherwise the result becomes 'inf'
    {
        return 1;
    }
}

void Tracer::OutputStatisticsPO(CalcTimer &timer, long long orNumber, const string &path)
{
    string startTime = ctime(&m_startTime);
    string totalTime = timer.Elapsed();
    time_t end = timer.Stop();
    string endTime = ctime(&end);

    m_summary += "\nStart of calculation = " + startTime
                 + "End of calculation   = " + endTime
                 + "\nTotal time of calculation = " + totalTime
                 + "\nTotal number of body orientation = " + to_string(orNumber);

    double passedEnergy = (m_handler->m_outputEnergy/m_incomingEnergy)*100;

    m_summary += "\nTotal incoming energy = " + to_string(m_incomingEnergy)
                 + "\nTotal outcoming energy (GO) = " + to_string(m_handler->m_outputEnergy)
                 + " (S/4 = " + to_string(m_particle->Area()/4)
                 + ")\nEnergy passed = " + to_string(passedEnergy) + '%';

    //	if (isNanOccured)
    //	{
    //		m_summary += "\n\nWARNING! NAN values occured. See 'log.txt'";
    //	}

    ofstream out(path + "\\out.dat", ios::out);

    out << m_summary;
    out.close();

    cout << m_summary;
}

void Tracer::SetIsOutputGroups(bool value)
{
    isOutputGroups = value;
}

//REF: объединить с предыдущим
void Tracer::TraceRandomPO2(int betaNumber, int gammaNumber, const ScatteringRange &bsCone,
                              const Tracks &tracks, double wave)
{
//	m_wavelength = wave;
//	CalcTimer timer;
//	long long count = 0;

//	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
//	ofstream outFile(m_resultDirName, ios::out);

//	vector<Beam> outBeams;
//	double beta, gamma;

//	double betaNorm = m_symmetry.beta/betaNumber;
//	double gammaNorm = m_symmetry.gamma/gammaNumber;

//	int halfGammaCount = gammaNumber/2;
//	int normIndex = gammaNorm/m_symmetry.gamma;

//	timer.Start();

//	for (int i = 0; i <= betaNumber; ++i)
//	{
//		beta = i*betaNorm;

//		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
//		{
//			gamma = j*gammaNorm;

//			for (size_t groupID = 0; groupID < tracks.size(); ++groupID)
//			{
//				bool ok = m_scattering->ScatterLight(beta, gamma, tracks[groupID].tracks, outBeams);

//				CleanJ(tracks.size(), bsCone);
//				HandleBeamsPO2(outBeams, bsCone, groupID);

//				outBeams.clear();
//				AddResultToMatrix(M, bsCone, normIndex);
//			}
//		}

//		WriteConusMatrices(outFile, M, bsCone);

//		OutputProgress(betaNumber, count, timer);
//		++count;
//	}

//	outFile.close();
}

void Tracer::OutputStartTime(CalcTimer &timer)
{
    m_startTime = timer.Start();
    cout << "Started at " << ctime(&m_startTime) << endl;
    cout << "Progress\tOrientation\tTime passed\tOutput file"<< endl;
}

void Tracer::SetHandler(Handler *handler)
{
    m_handler = handler;
    m_handler->SetScattering(m_scattering);
}

void Tracer::HandleBeamsPO2(vector<Beam> &outBeams, const ScatteringRange &bsCone, int groupID)
{
//	for (unsigned int i = 0; i < outBeams.size(); ++i)
//	{
//		Beam &beam = outBeams.at(i);

//		beam.RotateSpherical(-m_incidentLight.direction,
//							 m_incidentLight.polarizationBasis);

//		Point3f center = beam.Center();
//		Point3d center_d = Point3d(center);
//		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

//		Point3f T = CrossProduct(beam.polarizationBasis, beam.direction);
//		T = T/Length(T); // basis of beam

//		for (int i = 0; i <= bsCone.phiCount; ++i)
//		{
//			double f = i*bsCone.dPhi;
//			double sinF = sin(f), cosF = cos(f);

//			for (int j = 0; j <= bsCone.thetaCount; ++j)
//			{
//				double t = j*bsCone.dTheta;
//				double sinT = sin(t);

//				Point3d vr(sinT*cosF, sinT*sinF, cos(t));
//				Point3d vf = (j == 0) ? -m_incidentLight.polarizationBasis
//									  : Point3d(-sinF, cosF ,0);
//				matrixC Jn_rot(2, 2);
//				CalcJnRot(beam, T, vf, vr, Jn_rot);

//				complex fn(0, 0);
//				fn = beam.DiffractionIncline(vr, m_wavelength);

//				double dp = DotProductD(vr, center_d);
//				complex tmp = exp_im(M_2PI*(lng_proj0-dp)/m_wavelength);
//				matrixC fn_jn = beam.J * tmp;

//				matrixC c = fn*Jn_rot*fn_jn;
//				J[groupID].insert(i, j, c);
//			}
//		}
//	}
}
