#include "basics.h"

vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &ttau, const int expansion);

vector<vector<double> > RcPDF(function<double(double)> Gamma, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

vector<vector<double> > Nbar(function<double(double)> Gamma, const double x1, const double x2, const vector<vector<double> > &Ft, const vector<vector<double> > &taut, const vector<vector<double> > &at);

vector<double> findtrange(function<double(double)> Gamma, const double Nbarmin, const int Nb, const double tfrac, const double Fmin, const int expansion);

vector<bubble> nucleate(function<double(double)> Gamma, const double x1, const double x2, const int Nn, const vector<vector<double> > &taut, const vector<vector<double> > &at, rgen &mt);

double findtauc(const double x1, const double x2, const double taun, const vector<double> &xh, const vector<double> &xc, const vector<bubble> &bubbles, const int jb, const double taumax);

vector<complex<double> > PEEstep(int jt, vector<complex<double> > &udu, vector<vector<complex<double> > > &Tt, vector<vector<double> > &at, vector<vector<double> > &Ht, double k);

vector<complex<double> > TTprojection6(const vector<complex<double> > &X, const vector<double> k);
