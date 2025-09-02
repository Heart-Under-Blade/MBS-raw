#include <iostream>
#include <assert.h>
#include <float.h>
#include <string>

#include "HexagonalAggregate.h"
#include "ConcaveHexagonal.h"
#include "CertainAggregate.h"
#include "Bullet.h"
#include "BulletRosette.h"
#include "global.h"
#include "TracerGO.h"
#include "TracerPOTotal.h"
#include "ArgPP.h"
#include "Tracks.h"
#include "HandlerPOTotal.h"
#include "HandlerTotalGO.h"
#include "HandlerTracksGO.h"
#include "Droxtal.h"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;

enum class ParticleType : int
{
    Hexagonal = 1,
    Bullet = 2,
    BulletRosette = 3,
    Droxtal = 4,
    ConcaveHexagonal = 10,
    TiltedHexagonal = 11,
    HexagonalAggregate = 12,
    CertainAggregate = 999
};

Tracks trackGroups;

void SetArgRules(ArgPP &parser)
{
    int zero = 0;
    parser.AddRule("p", '+'); // particle (type, size, ...)
    parser.AddRule("ri", 2); // refractive index (Re and Im parts)
    parser.AddRule("n", 1); // number of internal reflection
    parser.AddRule("pf", zero, true); // particle (filename)
    parser.AddRule("rs", 1, true, "pf"); // resize particle (new size)
    parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
    parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
    parser.AddRule("go", 0, true); // geometrical optics method
    parser.AddRule("po", 0, true); // phisical optics method
    parser.AddRule("w", 1, true); // wavelength
    parser.AddRule("b", 2, true, "po"); // beta range (begin, end)
    parser.AddRule("g", 2, true, "po"); // gamma range (begin, end)
    parser.AddRule("grid", '+', true); /* backscattering grid:
 (radius, Nphi, Ntheta) when 3 parameters
 (theta1, theta2, Nphi, Ntheta) when 4 parameters*/
    parser.AddRule("point", zero, true, "po"); // calculate only backscatter point
    parser.AddRule("tr", 1, true); // file with trajectories
    parser.AddRule("all", 0, true); // calculate all trajectories
    parser.AddRule("abs", zero, true, "w"); // accounting of absorbtion
    parser.AddRule("close", 0, true); // closing of program after calculation
    parser.AddRule("o", 1, true); // output folder name
    parser.AddRule("gr", zero, true);
    parser.AddRule("filter", 1, true); // scattering angle filter
    parser.AddRule("shadow", zero, true);
    parser.AddRule("incoh", zero, true);
    parser.AddRule("shadow_off", zero, true);
    parser.AddRule("forced_nonconvex", zero, true);
    parser.AddRule("forced_convex", zero, true);
    parser.AddRule("r", 1, true); // restriction ratio for small beams when intersection (100 by default)
    parser.AddRule("log", 1, true); // time of writing progress (in seconds)
}

ScatteringRange SetConus(ArgPP &parser)
{
    if (parser.GetArgNumber("grid") == 3)
    {
        double radius = parser.GetDoubleValue("grid", 0);
        int nAz = parser.GetDoubleValue("grid", 1);
        int nZen = parser.GetDoubleValue("grid", 2);
        return ScatteringRange(M_PI - DegToRad(radius), M_PI, nAz, nZen);
    }
    else if (parser.GetArgNumber("grid") == 4)
    {
        double zenStart = parser.GetDoubleValue("grid", 0);
        double zenEnd = parser.GetDoubleValue("grid", 1);

        if (zenStart < zenEnd)
        {
            int nAz = parser.GetDoubleValue("grid", 2);
            int nZen = parser.GetDoubleValue("grid", 3);
            return ScatteringRange(DegToRad(zenStart), DegToRad(zenEnd), nAz, nZen);
        }
        else
        {
            std::cerr << "ERROR!!!!!!! In \"grid\" arg 1 must be less than arg 2";
            exit(1);
        }
    }
    else
    {
        std::cerr << "ERROR!!!!!!! Wrong \"grid\" argument number";
        exit(1);
    }
}

AngleRange GetRange(const ArgPP &parser, const std::string &key,
                    Particle *particle)
{
    int number;
    double min, max;

    if (key == "b")
    {
        number = parser.GetIntValue("random", 0);

        if (parser.IsCatched("b"))
        {
            min = DegToRad(parser.GetDoubleValue(key, 0));
            max = DegToRad(parser.GetDoubleValue(key, 1));
        }
        else
        {
            min = 0;
            max = particle->GetSymmetry().beta;
        }
    }
    else if (key == "g")
    {
        number = parser.GetIntValue("random", 1);

        if (parser.IsCatched("g"))
        {
            min = DegToRad(parser.GetDoubleValue(key, 0));
            max = DegToRad(parser.GetDoubleValue(key, 1));
        }
        else
        {
            min = 0;
            max = particle->GetSymmetry().gamma;
        }
    }
    else
    {
        cerr << "Error! " << __FUNCTION__;
        throw std::exception();
    }

    return AngleRange(min, max, number);
}

int main(int argc, const char* argv[])
{
    RenameConsole("MBS");

    std::string additionalSummary;

    if (argc <= 1) // no arguments
    {
        cout << "No arguments. Press any key to exit..." << endl;
        getchar();
        return 1;
    }

    ArgPP args;
    SetArgRules(args);
    args.Parse(argc, argv);

    bool isAbs = args.IsCatched("abs");

    double re = args.GetDoubleValue("ri", 0);
    double im = args.GetDoubleValue("ri", 1);
    complex refrIndex = complex(re, im);

    if (!isAbs)
    {
        refrIndex = complex(re, 0);
    }

    // TODO: AggregateBuilder

    Particle *particle = nullptr;


    additionalSummary += "Particle: ";

    if (args.IsCatched("pf"))
    {
        std::string filename = args.GetStringValue("p");
        particle = new Particle();
        particle->SetFromFile(filename);
        particle->SetRefractiveIndex(complex(refrIndex));

        double origDMax = particle->MaximalDimention();
        additionalSummary += "from file: " + filename + "\n";
        additionalSummary += "\tOriginal Dmax: " + std::to_string(origDMax);

        if (args.IsCatched("rs"))
        {
//            std::string dimension = args.GetStringValue("resize", 0);
            double newSize = args.GetDoubleValue("rs", 0);

//            if (dimension == "d")
//            {
                particle->Resize(newSize);
//            }
//            else if (dimension == "v")
//            {
//                double oldV = particle->Volume();
//                double ratio = newSize/oldV;
//                ratio = pow(ratio, 1.0/3);
//                particle->Scale(ratio);
//            }
        }

        double newDMax = particle->MaximalDimention();
        additionalSummary += ", new Dmax: " + std::to_string(newDMax)
                + ", resize factor: " + std::to_string(newDMax/origDMax) + '\n';
    }
    else
    {
        ParticleType type = (ParticleType)args.GetIntValue("p", 0);
        double height;
        double diameter;

        double sup;
        int num;

        switch (type)
        {
        case ParticleType::Hexagonal:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            particle = new Hexagonal(refrIndex, diameter, height);
            break;
        case ParticleType::Bullet:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = (diameter*sqrt(3)*tan(DegToRad(62)))/4;
            particle = new Bullet(refrIndex, diameter, height, sup);
            break;
        case ParticleType::BulletRosette:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = (args.GetArgNumber("p") == 4)
                    ? args.GetDoubleValue("p", 3)
                    : diameter*sqrt(3)*tan(DegToRad(62))/4;
            particle = new BulletRosette(refrIndex, diameter, height, sup);
            break;
        case ParticleType::Droxtal:
            sup = args.GetDoubleValue("p", 3);
            particle = new Droxtal(refrIndex, DegToRad(32.35), DegToRad(71.81), sup);
            break;
//		case ParticleType::TiltedHexagonal:
//			sup = parser.argToValue<double>(vec[3]);
//			particle = new TiltedHexagonal(r, hh, ri, sup);
//			break;
        case ParticleType::ConcaveHexagonal:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = args.GetDoubleValue("p", 3);
            particle = new ConcaveHexagonal(refrIndex, diameter, height, sup);
            break;
        case ParticleType::HexagonalAggregate:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            num = args.GetIntValue("p", 3);
            particle = new HexagonalAggregate(refrIndex, diameter, height, num);
            break;
        case ParticleType::CertainAggregate:
            sup = args.GetDoubleValue("p", 1);
            particle = new CertainAggregate(refrIndex, sup);
            break;
        default:
            assert(false && "ERROR! Incorrect type of particle.");
            break;
        }
    }

    additionalSummary += "\tRefractive index: " + to_string(re);

    if (fabs(im) > FLT_EPSILON)
    {
        additionalSummary +=  " + i\n";
    }

    additionalSummary += "\n";

    particle->Output("particle_for_check.dat");
    additionalSummary += "\tArea:" + to_string(particle->Area()) + "\n\n";

    int reflNum = args.GetDoubleValue("n");
    additionalSummary += "Number of secondary reflections: " + to_string(reflNum) + "\n";

    string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
                                           : "M";
    size_t pos = 0;

    while ((pos = dirName.find('%', pos)) != string::npos)
    {
        size_t start = pos;
        ++pos;
        string key;

        while (dirName[pos] != '_' && pos < dirName.size())
        {
            key += dirName[pos++];
        }

        if (key.size())
        {
            string val = args.GetStringValue(key);
            dirName.replace(start, pos-start, val);
            pos = start + val.size();
        }
    }

    bool isOutputGroups = args.IsCatched("gr");
    double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
    additionalSummary += "Wavelength (um): " + to_string(wave) + "\n";

    if (args.IsCatched("forced_nonconvex"))
    {
        particle->isConcave = true;
    }

    if (args.IsCatched("forced_convex"))
    {
        particle->isConcave = false;
    }

    if (args.IsCatched("tr"))
    {
        string trackFileName = args.GetStringValue("tr");
        trackGroups.ImportTracks(particle->nFacets, trackFileName);
        trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
        trackGroups.shouldOutputGroups = args.IsCatched("gr");
    }
    else
    {

    }

    int nTheta = args.GetDoubleValue("grid", 2);

    additionalSummary += "Method: ";

    if (args.IsCatched("po"))
    {
        additionalSummary += "Physical optics";

        if (args.IsCatched("fixed"))
        {
            additionalSummary += ", fixed orientation\n\n";

            ScatteringRange bsCone = SetConus(args);

            double beta  = args.GetDoubleValue("fixed", 0);
            double gamma = args.GetDoubleValue("fixed", 1);

            TracerPO tracer(particle, reflNum, dirName);
            tracer.m_summary = additionalSummary;

            HandlerPO *handler = new HandlerPO(particle, &tracer.m_incidentLight,
                                               nTheta, wave);
            if (args.IsCatched("r"))
            {
                tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
            }
            handler->isCoh = !args.IsCatched("incoh");
            handler->SetScatteringSphere(bsCone);
            handler->SetTracks(&trackGroups);
            handler->SetAbsorptionAccounting(isAbs);

            if (args.IsCatched("log"))
            {
                tracer.m_logTime = args.GetIntValue("log");
            }

            tracer.shadowOff = args.IsCatched("shadow_off");
            tracer.SetIsOutputGroups(isOutputGroups);
            tracer.SetHandler(handler);
            tracer.TraceFixed(beta, gamma);
        }
        else if (args.IsCatched("random"))
        {
            additionalSummary += ", random orientation\n\n";
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);

            HandlerPO *handler;

            if (args.IsCatched("point"))
            {
//                 TracerBackScatterPoint tracer(particle, reflNum, dirName);
//                 tracer.m_scattering->m_wave = wave;
//                 if (args.IsCatched("r"))
//                 {
//                     tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
//                 }
//                 cout << additionalSummary;
//                 tracer.m_summary = additionalSummary;

//                 handler = new HandlerBackScatterPoint(particle, &tracer.m_incidentLight,
//                                                       nTheta, wave);
//                 double normIndex = gamma.step/gamma.norm;
// //                handler->isCoh = !args.IsCatched("incoh");
//                 handler->SetNormIndex(normIndex);
//                 handler->SetTracks(&trackGroups);
//                 trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
//                 handler->SetAbsorptionAccounting(isAbs);

//                 if (args.IsCatched("filter"))
//                 {
//                     handler->SetBackScatteringConus(DegToRad(args.GetDoubleValue("filter")));
//                 }

//                 if (args.IsCatched("log"))
//                 {
//                     tracer.m_logTime = args.GetIntValue("log");
//                 }
//                 tracer.SetIsOutputGroups(isOutputGroups);
//                 tracer.SetHandler(handler);
//                 tracer.TraceRandom(beta, gamma);
            }
            else
            {
                TracerPO *tracer;
                ScatteringRange conus = SetConus(args);

                if (args.IsCatched("all"))
                {
                    tracer = new TracerPOTotal(particle, reflNum, dirName);
                    tracer->m_scattering->m_wave = wave;
                    if (args.IsCatched("r"))
                    {
                        tracer->m_scattering->restriction = args.GetDoubleValue("r", 0);
                    }
                    tracer->shadowOff = args.IsCatched("shadow_off");

                    trackGroups.push_back(TrackGroup());
                    handler = new HandlerPOTotal(particle, &tracer->m_incidentLight,
                                                 nTheta, wave);
                    handler->normIndexGamma = gamma.step/gamma.norm;

                    if (args.IsCatched("filter"))
                    {
                        handler->SetBackScatteringConus(DegToRad(args.GetDoubleValue("filter")));
                    }
                }
                else
                {
                    tracer = new TracerPO(particle, reflNum, dirName);
                    tracer->m_scattering->m_wave = wave;
                    tracer->shadowOff = args.IsCatched("shadow_off");
                    if (args.IsCatched("r"))
                    {
                        tracer->m_scattering->restriction = args.GetDoubleValue("r", 0);
                    }
                    handler = new HandlerPO(particle, &tracer->m_incidentLight,
                                            nTheta, wave);

                    if (args.IsCatched("filter"))
                    {
                        handler->SetBackScatteringConus(DegToRad(args.GetDoubleValue("filter")));
                    }
                }

                cout << additionalSummary;
                tracer->m_summary = additionalSummary;

                handler->isCoh = !args.IsCatched("incoh");
                handler->SetScatteringSphere(conus);
                handler->SetTracks(&trackGroups);
                handler->SetAbsorptionAccounting(isAbs);

                if (args.IsCatched("log"))
                {
                    tracer->m_logTime = args.GetIntValue("log");
                }
                tracer->SetIsOutputGroups(isOutputGroups);
                tracer->SetHandler(handler);
                tracer->TraceRandom(beta, gamma);
            }

            delete handler;
        }
        else
        {
            cout << endl << "error";
        }
    }
    else // go
    {
        additionalSummary += "Geometrical optics";

        TracerGO tracer(particle, reflNum, dirName);
        tracer.m_scattering->m_wave = wave;
        if (args.IsCatched("r"))
        {
            tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
        }
        tracer.m_summary = additionalSummary;
        tracer.SetIsOutputGroups(isOutputGroups);

        HandlerGO *handler;

        if (args.IsCatched("tr"))
        {
            handler = new HandlerTracksGO(particle, &tracer.m_incidentLight, nTheta, wave);
            handler->SetTracks(&trackGroups);
        }
        else
        {
            handler = new HandlerTotalGO(particle, &tracer.m_incidentLight, nTheta, wave);
        }

        tracer.m_summary = additionalSummary;
        ScatteringRange grid = SetConus(args);
        handler->isCoh = !args.IsCatched("incoh");
        handler->SetScatteringSphere(grid);
        handler->SetAbsorptionAccounting(isAbs);
        tracer.SetHandler(handler);

        if (args.IsCatched("fixed"))
        {
            additionalSummary += ", fixed orientation, ";
            double beta  = args.GetDoubleValue("fixed", 0);
            double gamma = args.GetDoubleValue("fixed", 1);
            additionalSummary += "zenith " + to_string(beta) + "°, azimuth " + to_string(gamma) + "°" + "\n\n";
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            tracer.TraceFixed(beta, gamma);
        }
        else if (args.IsCatched("random"))
        {
            additionalSummary += ", random orientation, ";
            AngleRange beta = GetRange(args, "b", particle);
            AngleRange gamma = GetRange(args, "g", particle);
            additionalSummary += "grid: " + to_string(beta.number) + "x" + to_string(gamma.number) + "\n\n";
            cout << additionalSummary;
            tracer.m_summary = additionalSummary;
            tracer.TraceRandom(beta, gamma);
        }

        delete handler;
    }

    cout << endl << "done";

    if (!args.IsCatched("close"))
    {
        getchar();
    }

    return 0;
}
