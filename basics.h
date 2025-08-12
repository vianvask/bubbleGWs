#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>

using namespace std;

typedef mt19937_64 rgen;

const double PI = 3.141592653589793238463;
const double beta = 1.0;
const complex<double> I(0.0, 1.0);
const complex<double> zero(0.0, 0.0);

class bubble {
    public:
        vector<double> x;
        double tau;
        double a;
};

string to_string_prec(const double a, const int n);

int delta(int i, int j);
int sgn(double a);
int imod(int i, int n);

double radius(double tau, double taun);
double randomreal(double x1, double x2, rgen &mt);
vector<vector<double> > sphereN(int Ns);
double inner(const vector<double> &X, const vector<double> &Y);
double distance(const vector<double> &X, const vector<double> &Y);
double distance(const vector<double> &X, const vector<double> &Y, double x1, double x2);
vector<double> periodic(double x1, double x2, vector<double> &x);
double interpolate(double x, vector<vector<double> > &y);
double findrootG(double y, double dx, vector<vector<double> > &list);
void writeToFile(vector<vector<double> > &matrix, const string &filename);
