#include "basics.h"

vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &atau, vector<vector<double> > &ttau);

vector<vector<double> > Nbar(function<double(double)> Gamma, const double x1, const double x2, const vector<vector<double> > &Ft, const vector<vector<double> > &taut, const vector<vector<double> > &at);

vector<double> findtrange(function<double(double)> Gamma, const double Nbarmin, const int Nb, const double tfrac);

vector<bubble> nucleate(const double x1, const double x2, const int Nn, const vector<vector<double> > &taut, const vector<vector<double> > &Nb, rgen &mt);

double findtauc(const double x1, const double x2, const double taun, const vector<double> &xh, const vector<double> &xc, const vector<bubble> &bubbles, const int jb, const double taumax);

vector<complex<double> > TTprojection6(const vector<complex<double> > &X, const vector<double> k);
