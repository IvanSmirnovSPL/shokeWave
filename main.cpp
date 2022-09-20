#include <iostream>
#include <cmath>

struct problemParms
{
    problemParms() = default;
    problemParms(double Mach, double beta): _Mach(Mach), _beta(beta){}
    double _Mach;
    double _beta;
};

struct waveParams
{
    waveParams(double pressure, double temperature, double rho):
    _rho(rho), _pressure(pressure), _tempreture(temperature){}
    waveParams(double pressure, double temperature):
    _pressure(pressure), _tempreture(temperature){_rho = M * _pressure / R / temperature;}
    waveParams() = default;
    double _rho;
    double _pressure;
    double _tempreture;
    static constexpr double R = 8.31446261815324; // universal gas const
    static constexpr double M = 29.e-3; //simple air
};

std::ostream& operator<<(std::ostream& out, const waveParams& P)
{
    out<<"Pressure: "<<P._pressure<<std::endl<<"Tempreture: "<<P._tempreture<<std::endl<<"Density: "<<P._rho<<std::endl;
    return out;
}

double methodNewton(double x_old, double pressision, double(*f)(double), double(*df)(double))
{  
    double x_new = x_old;
    bool logic = true;
    do
    {
        x_new = x_old - (*f)(x_old) / (*df)(x_old);
        logic = (fabs(x_old - x_new) > pressision);
        x_old = x_new;

    } while (logic);
    return x_new;
}

double f(double x)
{
    return x * x - 4 * x + 4;
}

double df(double x)
{
    return 2 * x  - 4;
}


int main()
{
    problemParms params(5., 15.);
    waveParams P1 = waveParams(1.01325e5, 300);
    double rez = methodNewton(1., 1.e-5, &f, &df);
    
    return 0;
}
