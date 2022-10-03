#include <iostream>
#include <cmath>

struct problemParms
{
    problemParms() = default;
    problemParms(double Mach, double beta): _Mach(Mach){_beta = beta * M_PI / 180.;}
    double _Mach;
    double _beta;
};

struct waveParams
{
    waveParams(const problemParms& PARAMS, double pressure, double temperature, double rho):
    _PARAMS(PARAMS), _rho(rho), _pressure(pressure), _tempreture(temperature){}
    waveParams(double pressure, double temperature, double rho):
    _rho(rho), _pressure(pressure), _tempreture(temperature){}
    waveParams(const problemParms& PARAMS, double pressure, double temperature):
    _PARAMS(PARAMS), _pressure(pressure), _tempreture(temperature){_rho = M * _pressure / R / temperature;
                                                  _a = std::sqrt(gamma * _pressure / _rho);
                                                  _v = PARAMS._Mach * _a;}
    waveParams() = default;

    waveParams afterWave(double theta)
    {
        using namespace std;
        _v_n = _v * sin(theta);
        _u_n = 1 - (sin(theta - _PARAMS._beta) * cos(theta)) /
                     (sin(theta) * cos(theta - _PARAMS._beta));
        double rho = _rho / (1 - _u_n);
        double pressure = _pressure + _rho * _u_n * _v_n * _v_n;
        double temperature = M * pressure / _rho / R;
        return waveParams(pressure, temperature, rho);
        

    }
    double _rho;
    double _pressure;
    double _tempreture;
    double _a;
    double _v;
    double _v_n, _u_n;
    problemParms _PARAMS;
    static constexpr double R = 8.31446261815324; // universal gas const
    static constexpr double M = 29.e-3; //simple air
    static constexpr double gamma = 1.4; 
};

std::ostream& operator<<(std::ostream& out, const waveParams& P)
{
    out<<"Pressure: "<<P._pressure<<std::endl<<"Tempreture: "<<P._tempreture<<std::endl<<"Density: "<<P._rho<<std::endl;
    return out;
}

double methodNewton(double x_old, double pressision, const waveParams& wP, const problemParms& pP,
                    double(*f)(double, const waveParams&, const problemParms&),
                    double(*df)(double, const waveParams&, const problemParms&))
{  
    double x_new = x_old;
    bool logic = true;
    do
    {
        x_new = x_old - (*f)(x_old, wP, pP) / (*df)(x_old, wP, pP);
        logic = (fabs(x_old - x_new) > pressision);
        x_old = x_new;

    } while (logic);
    return x_new;
}

double f(double theta, const waveParams& wP, const problemParms& pP)
{
    using namespace std;
    return sin(theta) * sin(theta) - 0.5 * (wP.gamma + 1) * sin(theta) * sin(pP._beta) / cos(theta - pP._beta)
                                                                             - (1 / pP._Mach) * (1 / pP._Mach);
}

double df(double theta, const waveParams& wP, const problemParms& pP)
{
    using namespace std;
    return 2 * sin(theta) * cos(theta) - 0.5 * (wP.gamma + 1) * cos(theta) * sin(pP._beta) / cos(theta - pP._beta) -
    0.5 * (wP.gamma + 1) * sin(theta - pP._beta) * sin(theta) * sin(pP._beta) / (cos(theta - pP._beta) * cos(theta - pP._beta));
}


int main()
{
    problemParms params(5., 15.);
    waveParams P1 = waveParams(params, 1.01325e5, 300);
    double theta = methodNewton(2 * params._beta + M_PI / 4, 1.e-5, P1, params, &f, &df);
    std::cout<<"theta: "<<theta * 180. / M_PI<<std::endl;
    std::cout<< P1.afterWave(theta)<<std::endl;
    return 0;
}
