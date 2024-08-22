#include "basics.h"

string to_string_prec(const double a, const int n) {
    ostringstream out;
    out.precision(n);
    out << fixed << a;
    return out.str();
}

int delta(int i, int j) {
    return (i==j);
}

int sgn(const double a) {
    return (0.0 < a) - (a < 0.0);
}

int imod(int i, int n) {
    return (i % n + n) % n;
}

// bubble radius as a function of conformal time
double radius(double tau, double taun) {
    double R = -1.0;
    if (tau >= taun)
        R = tau-taun;
    return R;
}

// random real number in the range (x_1,x_2)
double randomreal(double x1, double x2, rgen &mt) {
    long double r01 = mt()/(1.0*mt.max());
    return (x1 + (x2-x1)*r01);
}

// generates Ns semi-uniformly distributed points on the surface of an unit sphere
vector<vector<double> > sphereN(int Ns) {
    double theta, phi, phiold = 0.0, hk;
    vector<double> xh(3);
    vector<vector<double> > xhcoll(Ns, vector<double> (3));
    for (int jp = 0; jp < Ns; jp++) {
        hk = -1.0+2.0*jp/(Ns-1.0);
        theta = acos(hk);
        if (jp>0 && jp<Ns-1) {
            phi = fmod(phiold + 3.809/sqrt(1.0*Ns)/sqrt(1.0-hk*hk), 2.0*PI);
        } else {
            phi = 0.0;
        }
        phiold = phi;
        xh[0] = cos(phi)*sin(theta);
        xh[1] = sin(phi)*sin(theta);
        xh[2] = cos(theta);
        xhcoll[jp] = xh;
    }
    return xhcoll;
}

// distance between two points
double distance(const vector<double> &X, const vector<double> &Y) {
    double d2 = 0.0;
    for (int j = 0; j < X.size(); j++) {
        d2 += (X[j]-Y[j])*(X[j]-Y[j]);
    }
    return sqrt(d2);
}

// projection of vector x inside a cube
vector<double> periodic(double x1, double x2, vector<double> &x) {
    vector<double> y(3);
    for (int j=0; j<3; j++) {
        if (x[j] > 0) {
            y[j] = x[j] - x1;
        } else {
            y[j] = x[j] - x2;
        }
        
        y[j] = fmod(y[j],x2-x1);
        
        if (x[j] > 0) {
            y[j] = y[j] + x1;
        } else {
            y[j] = y[j] + x2;
        }
    }
    return y;
}

// linear interpolation
double interpolate(double x, vector<vector<double> > &y) {
    int n = y.size();
    if (x > y[n-1][0]) {
        cout << "Warning: the point lies above of the interpolation range." << endl;
        cout << x << "   " << y[n-1][0] << endl;
        return y[n-1][1];
    }
    if (x < y[0][0]) {
        cout << "Warning: the point lies below of the interpolation range." << endl;
        cout << x << "   " << y[0][0] << endl;
        return y[0][1];
    }
    double dx = y[1][0] - y[0][0];
    int jx = (int) ((x-y[0][0])/dx);
    if (jx < n-1) {
        return y[jx][1] + (y[jx+1][1] - y[jx][1])*(x - y[jx][0])/dx;
    }
    return y[n-1][1];
}
