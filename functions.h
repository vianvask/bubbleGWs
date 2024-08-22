#include "basics.h"

vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &atau);

vector<vector<double> > Nbar(function<double(double)> Gamma, const double x1, const double x2, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at);

vector<int> jtlist(vector<vector<double> > &Nk, int J, rgen &mt);

double findtauc(double x1, double x2, double taun, const vector<double> &xh, const vector<double> &xc, const vector<bubble> &bubbles, int jb, double taumax);

vector<complex<double> > TTprojection6(vector<complex<double> > &X, vector<double> k);
