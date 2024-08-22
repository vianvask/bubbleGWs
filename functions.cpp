#include "functions.h"

// evolution of the Universe on average, returns kmax
vector<double> averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &atau) {
    
    // initial state in vacuum dominance:
    double t = tmin, a = exp(tmin), tau = 1.0 - exp(-tmin);
    double rhoV = 1.0, rhoV0 = 1.0, rhoR = 0.0;
    
    double H = 1.0, F = 1.0, F0 = 1.0;
    vector<double> tmp(2);
    double Nt, kmax = 0.0, tkmax;
    
    for (int jt = 0; jt < jtmax; jt++) {
        tmp[0] = t;
        tmp[1] = F;
        Ft.push_back(tmp);
        tmp[1] = H;
        Ht.push_back(tmp);
        tmp[1] = a;
        at.push_back(tmp);
        tmp[1] = tau;
        taut.push_back(tmp);
        
        tmp[0] = tau;
        tmp[1] = a;
        atau.push_back(tmp);
        
        //calculate the false vacuum fraction
        Nt = 0.0;
        F0 = F;
        for (int j = 0; j < taut.size(); j++) {
            Nt += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*radius(tau,taut[j][1]), 3.0);
        }
        F = exp(-Nt);
        
        // compute the scale factor and conformal time
        H = sqrt(rhoV + rhoR);
        a += dt*H*a;
        tau += dt/a;
        
        if (a*H > kmax) {
            kmax = a*H;
            tkmax = t;
        }
        
        // update the energy densities
        rhoV = F;
        rhoR += -4.0*H*rhoR*dt - (F-F0);
        
        t += dt;
    }
    
    tmp[0] = kmax;
    tmp[1] = tkmax;
    
    return tmp;
}

// expected number of bubbles in sphere of radius 1/k
vector<vector<double> > Nbar(function<double(double)> Gamma, const double x1, const double x2, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    
    double t, tau, Np, N = 0.0;
    vector<double> tmp(3);
    vector<vector<double> > Nt;
    for (int jt = 0; jt < at.size(); jt++) {
        t = at[jt][0];
        tau = taut[jt][1];
        
        Np = N;
        N = 0.0;
        for (int j = 0; j < jt; j++) {
            N += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*(x2-x1),3.0);
        }
        
        tmp[0] = t;
        tmp[1] = N;
        tmp[2] = (N-Np)/dt;
        Nt.push_back(tmp);
    }
    
    return Nt;
}

// generate times tj for j<J
vector<int> jtlist(vector<vector<double> > &Nk, int J, rgen &mt) {
    const double dt = Nk[1][0] - Nk[0][0];
    vector<int> jtlist(J, Nk.size() - 1);
    int j = 0;
    for (int jt = 0; jt < Nk.size(); jt++) {
        if (dt*Nk[jt][2] > 1.0) {
            cout << "Warning: too long nucleation timestep, N = " << Nk[jt][1] << "." << endl;
        }
        if (dt*Nk[jt][2] > randomreal(0.0,1.0,mt)) {
            jtlist[j] = jt;
            j++;
        }
        if (j >= J) {
            return jtlist;
        }
    }
    cout << "Warning: nucleation finished at j = " << j << "." << endl;
    return jtlist;
}

// finds the conformal time when the bubble that nucleated at xc collided in the direction xh with some other bubble
double findtauc(double x1, double x2, double taun, const vector<double> &xh, const vector<double> &xc, const vector<bubble> &bubbles, int jb, double taumax) {
    
    double tauc = taumax, taun2, d, ct;
    vector<double> xc2(3), xc3(3);
    for (int j = 0; j < bubbles.size(); j++) {
        taun2 = bubbles[j].tau;
        xc2 = bubbles[j].x;
        if (j!=jb) {
            // account for periodic boundary considitions
            for (int j1 = -1; j1 < 2; j1++) {
                for (int j2 = -1; j2 < 2; j2++) {
                    for (int j3 = -1; j3 < 2; j3++) {
                        xc3[0] = xc2[0] + j1*(x2-x1);
                        xc3[1] = xc2[1] + j2*(x2-x1);
                        xc3[2] = xc2[2] + j3*(x2-x1);
                        d = distance(xc, xc3);
                        ct = (xh[0]*(xc3[0]-xc[0]) + xh[1]*(xc3[1]-xc[1]) + xh[2]*(xc3[2]-xc[2]))/d;
                        
                        // calculate the collision time using cosine rule
                        if (ct*d+(taun-taun2) > 0) {
                            tauc = min(tauc, taun + (pow(d,2.0)-pow(taun-taun2,2.0))/(2.0*(ct*d+(taun-taun2))));
                        }
                    }
                }
            }
        } else {
            // check self-collisions in periodic box
            for (int j1 = -1; j1 < 2; j1++) {
                for (int j2 = -1; j2 < 2; j2++) {
                    for (int j3 = -1; j3 < 2; j3++) {
                        if (!(j1==0 && j2==0 && j3==0)) {
                            xc3[0] = xc2[0] + j1*(x2-x1);
                            xc3[1] = xc2[1] + j2*(x2-x1);
                            xc3[2] = xc2[2] + j3*(x2-x1);
                            d = distance(xc, xc3);
                            ct = (xh[0]*(xc3[0]-xc[0]) + xh[1]*(xc3[1]-xc[1]) + xh[2]*(xc3[2]-xc[2]))/d;
                            
                            // calculate the collision time using cosine rule
                            if (ct*d+(taun-taun2) > 0) {
                                tauc = min(tauc, taun + (pow(d,2.0)-pow(taun-taun2,2.0))/(2.0*(ct*d+(taun-taun2))));
                            }
                        }
                    }
                }
            }
        }
    }
    return tauc;
}

vector<vector<complex<double> > > TTprojection(vector<vector<complex<double> > > &X, vector<double> k) {
    vector<vector<complex<double> > > Y(3, vector<complex<double> > (3, zero));
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int l=0; l<3; l++) {
                for (int m=0; m<3; m++) {
                    Y[i][j] = (1.0*delta(i,l)*delta(j,m) - 2.0*delta(i,l)*k[j]*k[m] + 0.5*delta(i,j)*k[l]*k[m] + 0.5*delta(l,m)*k[i]*k[j] - 0.5*delta(i,j)*delta(l,m) + 0.5*k[i]*k[j]*k[l]*k[m])*X[l][m];
                }
            }
        }
    }
    return Y;
}

vector<complex<double> > TTprojection6(vector<complex<double> > &X, vector<double> k) {
    vector<vector<complex<double> > > Y(3, vector<complex<double> > (3, zero));
    vector<complex<double> > Z(6, zero);
    
    Y[0][0] = X[0];
    Y[0][1] = X[1];
    Y[0][2] = X[2];
    Y[1][1] = X[3];
    Y[1][2] = X[4];
    Y[2][2] = X[5];
    
    Y[1][0] = Y[0][1];
    Y[2][0] = Y[0][2];
    Y[2][1] = Y[1][2];

    Y = TTprojection(Y, k);
    
    Z[0] = Y[0][0];
    Z[1] = Y[0][1];
    Z[2] = Y[0][2];
    Z[3] = Y[1][1];
    Z[4] = Y[1][2];
    Z[5] = Y[2][2];
    
    return Z;
}

