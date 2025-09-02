#pragma once

#include "geometry_lib.h"
#include "CalcTimer.h"
#include "Handler.h"

struct AngleRange
{
    double min;
    double max;
    int number;
    double norm;
    double step;

    AngleRange(double _min, double _max, int _number)
        : number(_number)
    {
        min = _min;
        max = _max;
        norm = max - min;
        step = norm/number;
    }
};

class Tracer
{
public:
    Tracer(Particle *particle, int nActs, const std::string &resultFileName);
    ~Tracer();

    // REF: delete?
    void TraceRandomPO2(int betaNumber, int gammaNumber, const ScatteringRange &bsCone,
                        const Tracks &tracks, double wave);

    void SetHandler(Handler *handler);

    void SetIsOutputGroups(bool value);// REF: заменить

    void CalcCsBeta(int betaNorm, double beta, const AngleRange &betaRange,
                    const AngleRange &gammaRange, double normIndex, double &cs_beta);

    void OutputStatisticsPO(CalcTimer &timer, long long orNumber, const std::string &path);
    double CalcNorm(long long orNum);

    bool shadowOff = false;
    Light m_incidentLight;
    std::string m_summary;
    int m_logTime = 5;
    Scattering *m_scattering;

protected:
    Handler *m_handler;
    Particle *m_particle;

    double m_incomingEnergy;
    double m_outcomingEnergy;

    std::string m_resultDirName;
    double m_wavelength;
    Symmetry m_symmetry;
    time_t m_startTime;

    // REF: заменить
    bool isOutputGroups = false;

    long long m_timeElapsed = 0;

protected:
    void OutputStartTime(CalcTimer &timer);
    void OutputProgress(int nOrientation, long long count,
                        int zenith, int azimuth, CalcTimer &timer, int nBeams);
    void OutputOrientationToLog(int i, int j, std::ostream &logfile);

private:
    void HandleBeamsPO2(std::vector<Beam> &outBeams, const ScatteringRange &bsCone, int groupID);
    void SetIncidentLight(Particle *particle);
};
