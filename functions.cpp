#include "functions.h"

// evolution of the Universe on average, returns (k_max, t_{k_max}) with expansion and (a_p, t_p) without
double averageevolution(function<double(double)> Gamma, const double tmin, const int jtmax, const double dt, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at, vector<vector<double> > &Ht, vector<vector<double> > &ttau, const int expansion) {
    
    // initial state in vacuum dominance:
    double t = tmin, a = exp(tmin), tau = 1.0 - exp(-tmin);
    double rhoV = 1.0, rhoV0 = 1.0, rhoR = 0.0;
    double H = 1.0, F = 1.0, F0 = 1.0;
    
    if (expansion == 0) {
        a = 1.0;
        tau = tmin;
        H = 0.0;
    }
    
    vector<double> tmp(2);
    double Nt, tp;
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
        tmp[1] = t;
        ttau.push_back(tmp);
        
        // compute the false vacuum fraction
        Nt = 0.0;
        F0 = F;
        for (int j = 0; j < taut.size(); j++) {
            Nt += 4.0*PI/3.0*dt*Gamma(taut[j][0])*pow(at[j][1]*radius(tau,taut[j][1]), 3.0);
        }
        F = exp(-Nt);
        
        // compute the scale factor and conformal time
        H = sqrt(rhoV + rhoR);
        if (expansion == 0) {
            H = 0.0;
        }
        a += dt*H*a;
        tau += dt/a;
        
        if (F > 0.3) {
            tp = t;
        }
        
        // update the energy densities
        rhoV = F;
        rhoR += -4.0*H*rhoR*dt - (F-F0);
        
        t += dt;
    }
    return tp;
}

// compute the distribuion of collision radii
vector<vector<double> > RcPDF(function<double(double)> Gamma, vector<vector<double> > &Ft, vector<vector<double> > &taut, vector<vector<double> > &at) {
    vector<vector<double> > pRc(50000, vector<double> (2,0.0));
    double dR = (taut.back()[1] - taut.front()[1])/(1.0*pRc.size());
    
    // compute the time derivative of F
    const double dt = at[1][0] - at[0][0];
    vector<vector<double> > DFtau(Ft.size(), vector<double> (2,0.0));
    DFtau[0][0] = taut[0][1];
    for (int jt = 1; jt < at.size(); jt++) {
        DFtau[jt][0] = taut[jt][1];
        DFtau[jt][1] = (Ft[jt][1] - Ft[jt-1][1])/dt;
    }
    vector<vector<double> > Ftau(at.size(), vector<double> (2,0.0));
    vector<vector<double> > atau(at.size(), vector<double> (2,0.0));
    for (int jt = 0; jt < at.size(); jt++) {
        Ftau[jt][0] = taut[jt][1];
        Ftau[jt][1] = Ft[jt][1];
        
        atau[jt][0] = taut[jt][1];
        atau[jt][1] = at[jt][1];
    }
    
    double Rc, tn, an, etan, p, norm = 0.0;
    for (int jR = 0; jR < pRc.size(); jR++) {
        Rc = jR*dR;
        
        // compute integral over time
        p = 0.0;
        for (int jt = 0; jt < at.size()-1; jt++) {
            tn = at[jt][0];
            an = at[jt][1];
            etan = taut[jt][1];
            
            if (etan + Rc < DFtau.back()[0]) {
                p += -pow(an,3.0)*Gamma(tn)*interpolate(etan + Rc, atau)*interpolate(etan + Rc, Ftau)*interpolate(etan + Rc, DFtau)*dt;
            }
        }
        
        pRc[jR][0] = Rc;
        pRc[jR][1] = p;
        
        norm += p*dR;
    }
    
    // normalize
    for (int jR = 0; jR < pRc.size(); jR++) {
        pRc[jR][1] = pRc[jR][1]/norm;
    }
    
    return pRc;
}


// expected number of bubbles nucleated in cube [x_1,x_2]^3
vector<vector<double> > Nbar(function<double(double)> Gamma, const double x1, const double x2, const vector<vector<double> > &Ft, const vector<vector<double> > &taut, const vector<vector<double> > &at) {
    
    const double dt = at[1][0] - at[0][0];
    
    double t, Np, N = 0.0;
    vector<double> tmp(3);
    vector<vector<double> > Nt;
    for (int jt = 0; jt < at.size(); jt++) {
        t = at[jt][0];
        
        Np = N;
        N = 0.0;
        for (int j = 0; j < jt; j++) {
            N += dt*Ft[j][1]*Gamma(at[j][0])*pow(at[j][1]*(x2-x1),3.0);
        }
        
        tmp[0] = t;
        tmp[1] = N;
        tmp[2] = (N-Np)/dt;
        Nt.push_back(tmp);
    }
    
    return Nt;
}

// finds the time range where the computation should be performed as well as the simulation boundaries
vector<double> findtrange(function<double(double)> Gamma, const double Nbarmin, const int Nb, const double tfrac, const double Fmin, const int expansion) {
    vector<double> trange(4);
    
    vector<vector<double> > Ft, taut, at, Ht, ttau;
    vector<double> tmp(2);
    
    int jtmax = 8000;
    double dt = 0.001;
        
    double tp = averageevolution(Gamma, -4.0, jtmax, dt, Ft, taut, at, Ht, ttau, expansion);
    
    vector<vector<double> > Tt(jtmax, vector<double> (2,0.0));
    for (int jt = 0; jt < jtmax; jt++) {
        Tt[jt][0] = Ft[jt][0];
        Tt[jt][1] = 1-Ft[jt][1];
    }

    vector<vector<double> > Nk = Nbar(Gamma, -0.5, 0.5, Ft, taut, at);
    int jtmaxnum = 0;
    for (int jt = 0; jt < jtmax; jt++) {
        if (isfinite(Nk[jt][1])) {
            jtmaxnum += 1;
        }
    }
    vector<vector<double> > Nt(jtmaxnum, vector<double> (2,0.0));
    for (int jt = 0; jt < jtmaxnum; jt++) {
        Nt[jt][0] = Nk[jt][0];
        Nt[jt][1] = Nk[jt][1];
    }
    
    double L = pow(Nb/Nt.back()[1], 1.0/3.0); // bar{N}(t=t_max) = N_b, fixes L
    
    trange[0] = findrootG(Nbarmin/pow(L,3.0), dt, Nt); // bar{N}(t=t_min) = barNtmin, fixes t_min
    trange[1] = tp + tfrac*L/pow(4.0*PI/3.0*Nb,1.0/3.0); // t_max = t_p + ftmax*<R>, fixes t_max
    trange[2] = L;
    trange[3] = findrootG(1.0-Fmin, dt, Tt); // bar{F}(t=t_max,nuc) = barFtmaxnuc, fixes t_max,nuc
    
    return trange;
}

// nucleates bubbles
vector<bubble> nucleate(function<double(double)> Gamma, const double x1, const double x2, const int Nn, const vector<vector<double> > &taut, const vector<vector<double> > &at, rgen &mt) {
    const double dt = taut[1][0] - taut[0][0];
    
    // generate list of nucleation sites
    vector<vector<double> > xlist(Nn, vector<double> (3, 0.0));
    for (int j = 0; j < Nn; j++) {
        for (int i = 0; i < 3; i++) {
            xlist[j][i] = randomreal(x1,x2,mt);
        }
    }
    
    bubble b0;
    vector<bubble> bubbles;
    double tn, taun, an;
    vector<double> xc(3, 0.0);
    
    // nucleate bubbles
    int Ntry = 0;
    const double dVcf = pow(x2-x1,3.0)/(1.0*Nn);
    bool flag, toohighprob = false;
    double dN;
    for (int jt = 0; jt < taut.size(); jt++) {
        tn = taut[jt][0];
        taun = taut[jt][1];
        an = at[jt][1];
        for (int j = 0; j < Nn; j++) {
            xc = xlist[j];
            
            dN = dVcf*dt*Gamma(tn)*pow(an,3.0);
            if (dN > 1.0) {
                toohighprob = true;
            }
            
            // try to nucleate a bubble
            if (dN > randomreal(0.0,1.0,mt)) {
                                
                // check if the bubble is inside another bubble
                flag = true;
                for (int jb = 0; jb < bubbles.size(); jb++) {
                    if (distance(xc, bubbles[jb].x, x1, x2) < radius(taun, bubbles[jb].tau)) {
                        flag = false;
                        jb = bubbles.size();
                    }
                }
                
                // nucleate the bubble only if it is not inside another bubble
                if (flag) {
                    b0.x = xc;
                    b0.tau = taun;
                    b0.a = an;
                    bubbles.push_back(b0);
                }
            }
        }
    }
    
    if (toohighprob) {
        cout << "Warning: too high nucleation probability." << endl;               
    }
    return bubbles;
}

// finds the conformal time when the bubble that nucleated at xc collided in the direction xh with some other bubble
double findtauc(const double x1, const double x2, const double taun, const vector<double> &xh, const vector<double> &xc, const vector<bubble> &bubbles, const int jb, const double taumax) {
    
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

// perturbed Einstein equation
complex<double> ddu(complex<double> u, complex<double> du, complex<double> T, double a, double H, double k) {
    return T/pow(a,2.0) - 3.0*H*du - pow(k/a,2.0)*u;
}

// take a step according to the perturbed Einstein equation using the 4th order Runge-Kutta method
vector<complex<double> > PEEstep(int jt, vector<complex<double> > &udu, vector<vector<complex<double> > > &Tt, vector<vector<double> > &at, vector<vector<double> > &Ht, double k) {
    
    complex<double> u = udu[0];
    complex<double> du = udu[1];
    
    double dt = at[jt+2][0] - at[jt][0];
    complex<double> T;
    double a, H;
    
    // t
    T = Tt[jt][1];
    a = at[jt][1];
    H = Ht[jt][1];
    complex<double> k11 = dt*du;
    complex<double> k12 = dt*ddu(u,du,T,a,H,k);
    
    // t + dt/2
    T = Tt[jt+1][1];
    a = at[jt+1][1];
    H = Ht[jt+1][1];
    complex<double> k21 = dt*(du+k12/2.0);
    complex<double> k22 = dt*ddu(u+k11/2.0,du+k12/2.0,T,a,H,k);
    complex<double> k31 = dt*(du+k22/2.0);
    complex<double> k32 = dt*ddu(u+k21/2.0,du+k22/2.0,T,a,H,k);

    // t + dt
    T = Tt[jt+2][1];
    a = at[jt+2][1];
    H = Ht[jt+2][1];
    complex<double> k41 = dt*(du+k32);
    complex<double> k42 = dt*ddu(u+k31,du+k32,T,a,H,k);

    vector<complex<double> > udu2(2,zero);
    udu2[0] = u + (k11+2.0*k21+2.0*k31+k41)/6.0;
    udu2[1] = du + (k12+2.0*k22+2.0*k32+k42)/6.0;
    
    return udu2;
}

// transvese traceless projection
vector<vector<complex<double> > > TTprojection(const vector<vector<complex<double> > > &X, const vector<double> khat) {
    vector<vector<complex<double> > > Y(3, vector<complex<double> > (3, zero));
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            // sum over l and m
            for (int l=0; l<3; l++) {
                for (int m=0; m<3; m++) {
                    Y[i][j] += (delta(i,l)*delta(j,m) - delta(i,j)*delta(l,m)/2.0 - delta(i,l)*khat[j]*khat[m] - delta(j,m)*khat[i]*khat[l] + delta(i,j)*khat[l]*khat[m]/2.0 + delta(l,m)*khat[i]*khat[j]/2.0 + khat[i]*khat[j]*khat[l]*khat[m]/2.0)*X[l][m];
                }
            }
        }
    }
    return Y;
}

// transvese traceless projection with input including only the 6 free components
vector<complex<double> > TTprojection6(const vector<complex<double> > &X, const vector<double> khat) {
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

    Y = TTprojection(Y, khat);
    
    Z[0] = Y[0][0];
    Z[1] = Y[0][1];
    Z[2] = Y[0][2];
    Z[3] = Y[1][1];
    Z[4] = Y[1][2];
    Z[5] = Y[2][2];
    
    return Z;
}

