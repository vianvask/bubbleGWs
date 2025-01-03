#include "functions.h"

int main (int argc, char *argv[]) {
    
    // input
    const double beta = atof(argv[1]);
    const double gammapbeta = atof(argv[2]);
    const int index = atoi(argv[3]);
    const int expansion = atoi(argv[4]);
    
    cout << setprecision(5) << fixed;
    clock_t time_req = clock(); // timing
    rgen mt(time(NULL)*(index+1)/beta); // random number generator

    // bubble nucleation rate, units chosen such that H0 = 1
    function<double(double)> Gamma = [beta, gammapbeta](double t) {
        return exp(beta*t - pow(gammapbeta*beta*t,2.0)/2.0);
    };
    
    int Nn = 10000; // #nucleation sites
    int Ns = 4000; // #points on the bubble surfaces
    int Nt = 8000; // #timesteps
    int Nk = 120; // #k values, k_max = Nk*k_min

    int J = 16; // bar{N}(t=t_p) = J, fixes L
    double barNtmin = 0.01; // bar{N}(t=t_min) = barNtmin, fixes t_min
    double ftmax = 4.0; // bar{F}(t=t_max) = t_p + ftmax*(t_p - t_1), fixes t_max
    double barFtmaxnuc = 0.001; // bar{F}(t=t_max,nuc) = barFtmaxnuc, fixes t_max,nuc

    // determine the time range and simulation volume
    vector<double> trange = findtrange(Gamma, barNtmin, J, ftmax, barFtmaxnuc, expansion);
    double tmin = trange[0];
    double tmax = trange[1];
    double tmaxnuc = min(trange[3],tmax);
    double dt = (tmax-tmin)/(1.0*Nt);
    double dtnuc = (tmaxnuc-tmin)/(1.0*Nt);
    double L = trange[2];
    double x1 = -L/2.0, x2 = L/2.0;
    
    // evolution of conformal time, scale factor, Hubble rate and expected number of bubbles, Nbar[jt][N,dN/dt]
    vector<vector<double> > Ft, taut, at, Ht, atau, ttau, Nb;
    averageevolution(Gamma, tmin, Nt, dtnuc, Ft, taut, at, Ht, atau, ttau, expansion);
    Nb = Nbar(Gamma, x1, x2, Ft, taut, at);
    
    // nucleate bubbles
    vector<bubble> bubbles;
    int Ntry = 0;
    while ((bubbles.size() < 8 || bubbles.size() > J) && Ntry < 100) {
        bubbles = nucleate(x1, x2, Nn, taut, Nb, mt);
        cout << bubbles.size() << endl;
        Ntry++;
    }
    cout << "#bubbles = " << bubbles.size() << endl;
    
    // recompute the background evolution with up to larger times with a larger dt
    Ft.clear(); taut.clear(); at.clear(); Ht.clear(); atau.clear(); ttau.clear();
    averageevolution(Gamma, tmin, Nt, dt, Ft, taut, at, Ht, atau, ttau, expansion);
    double taumax = taut[taut.size()-1][1];
    
    // generate a list of k values
    double kmin = 1.0/(2.0*L);
    double k = kmin;
    vector<double> klist;
    while (k < (Nk + 0.1)*kmin) {
        klist.push_back(k);
        k += kmin;
    }
    
    // generate x and k directions
    vector<vector<double> > xhat = sphereN(Ns);
    vector<vector<double> > khat(3, vector<double> (3, 0.0));
    int Nkhat = khat.size();
    for (int jd = 0; jd < Nkhat; jd++) {
        for (int j = 0; j < 3; j++) {
            khat[jd][j] = delta(jd, j);
        }
    }
    
    // initialize T_ij, u_ij and its time derivative: Nt timesteps, Nkhat k directions, Nk k values, Na approximations, 6 ij components
    int Na = 3;
    vector<vector<vector<vector<vector<complex<double> > > > > > T(Nt, vector<vector<vector<vector<complex<double> > > > > (Nkhat, vector<vector<vector<complex<double> > > > (Nk, vector<vector<complex<double> > > (Na, vector<complex<double> > (6, zero)))));
    vector<vector<vector<vector<vector<complex<double> > > > > > u = T, du = T;
    
    string filename = "B_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileB;
    outfileB.open(filename.c_str());
    
    double tau, taun, tauc, t, tc, a, ac, H, R, Rc, dV, kX, theta, phi;
    vector<double> xc(3, 0.0), xh(3, 0.0), X(3, 0.0), xixj(6, 0.0), F(Na, 0.0);

    // loop over all the bubbles
    for (int jb = 0; jb < bubbles.size(); jb++) {
        cout << jb << endl;
        
        xc = bubbles[jb].x;
        taun = bubbles[jb].tau;
        
        // loop over the points on the sphere
        for (int jp = 0; jp < Ns; jp++) {
            xh = xhat[jp];
            
            // convention as in TTprojection6
            xixj[0] = xh[0]*xh[0];
            xixj[1] = xh[0]*xh[1];
            xixj[2] = xh[0]*xh[2];
            xixj[3] = xh[1]*xh[1];
            xixj[4] = xh[1]*xh[2];
            xixj[5] = xh[2]*xh[2];
            
            // find when the collision happens in the direction xh
            tauc = findtauc(x1, x2, taun, xh, xc, bubbles, jb, 2.0*taumax);
            tc = interpolate(tauc, ttau);
            ac = interpolate(tauc, atau);
            Rc = radius(tauc, taun);
            
            // output the collision radii of first bubble
            if (jb < 1) {
                theta = acos(xh[2]);
                phi = sgn(xh[1]);
                if (sqrt(xh[0]*xh[0]+xh[1]*xh[1]) > 0) {
                    phi = sgn(xh[1])*acos(xh[0]/sqrt(xh[0]*xh[0]+xh[1]*xh[1]));
                }
                outfileB << theta << "   " << phi << "   " << Rc << "   " << ac << endl;
            }
            
            // compute the contribution to the stress-energy tensor
            for (int jt = 0; jt < Nt; jt++) {
                tau = taut[jt][1];
                if (tau > taun) {
                    R = radius(tau, taun);
                    dV = 4.0*PI/(3.0*Ns)*pow(R,3.0);
                    for (int j = 0; j < 3; j++) {
                        X[j] = xc[j] + R*xh[j];
                    }
                                        
                    // 0: envelope, 1: xi = 2, 2: xi = 3
                    if (R < Rc) {
                        F[0] = 1.0;
                        F[1] = 1.0;
                        F[2] = 1.0;
                    } else {
                        F[0] = 0.0;
                        F[1] = pow(Rc/R,3.0);
                        F[2] = pow(Rc/R,4.0);
                    }
                    
                    for (int jk = 0; jk < Nk; jk++) {
                        k = klist[jk];
                        for (int jd = 0; jd < Nkhat; jd++) {
                            kX = k*inner(khat[jd],X);
                            for (int ja = 0; ja < Na; ja++) {
                                for (int j6 = 0; j6 < 6; j6++) {
                                    T[jt][jd][jk][ja][j6] += dV*F[ja]*xixj[j6]*exp(-I*kX);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    outfileB.close();
    
    // solve the perturbed Einstein equation
    complex<double> u0, um1, du0, T0;
    for (int jt = 1; jt < Nt-1; jt++) {
        a = at[jt][1];
        H = Ht[jt][1];
        
        for (int jk = 0; jk < Nk; jk++) {
            k = klist[jk];
            for (int jd = 0; jd < Nkhat; jd++) {
                for (int ja = 0; ja < Na; ja++) {
                    for (int j6 = 0; j6 < 6; j6++) {
                        T0 = T[jt][jd][jk][ja][j6];
                        u0 = u[jt][jd][jk][ja][j6];
                        um1 = u[jt-1][jd][jk][ja][j6];
                        
                        u[jt+1][jd][jk][ja][j6] = (2.0*pow(dt/a,2.0)*(T0 - pow(k,2.0)*u0) + 4.0*u0 - 2.0*um1 + 3.0*dt*H*um1)/(2.0 + 3.0*dt*H);
                    }
                }
            }
        }
    }
    
    for (int jk = 0; jk < Nk; jk++) {
        for (int jd = 0; jd < Nkhat; jd++) {
            for (int ja = 0; ja < Na; ja++) {
                
                // compute the TT projection of u_ij(k)
                for (int jt = 0; jt < Nt; jt++) {
                    u[jt][jd][jk][ja] = TTprojection6(u[jt][jd][jk][ja], khat[jd]);
                }
                
                // compute the derivative of u_ij(k)
                for (int jt = 1; jt < Nt-1; jt++) {
                    for (int j6 = 0; j6 < 6; j6++) {
                        du[jt][jd][jk][ja][j6] = (u[jt+1][jd][jk][ja][j6] - u[jt-1][jd][jk][ja][j6])/(2.0*dt);
                    }
                }
            }
        }
    }
    Nt = Nt-1;
    
    filename = "OmegaGW_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmega;
    outfileOmega.open(filename.c_str());
    filename = "OmegaTotGW_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmegaTot;
    outfileOmegaTot.open(filename.c_str());
    
    const vector<double> zeroNa(Na, 0.0);
    const vector<double> j6coef {1.0, 2.0, 2.0, 1.0, 2.0, 1.0};
    vector<double> Omega(Na), OmegaTot(Na);
    double Theta;
    
    // output the GW spectrum and the total GW energy density as a function of time
    for (int jt = 0; jt < Nt; jt++) {
        t = at[jt][0];
        a = at[jt][1];
        H = Ht[jt][1];
        if (expansion > 0) {
            Theta = 4.0*PI/(1.0*Nkhat)*3.0*pow(1.0/H,2.0)/(16.0*pow(PI*L,3.0));
        } else {
            Theta = 4.0*PI/(1.0*Nkhat)*3.0/(16.0*pow(PI*L,3.0));
        }
        OmegaTot = zeroNa;
        for (int jk = 0; jk < Nk; jk++) {
            k = klist[jk];
            Omega = zeroNa;
            for (int ja = 0; ja < Na; ja++) {
                for (int jd = 0; jd < Nkhat; jd++) {
                    for (int j6 = 0; j6 < 6; j6++) {
                        du0 = du[jt][jd][jk][ja][j6];
                        Omega[ja] += pow(k,3.0)*Theta*j6coef[j6]*pow(abs(du0),2.0);
                    }
                    OmegaTot[ja] += kmin/(3.0*k)*Omega[ja];
                }
            }
            if (jt%20 == 0 || jt == Nt-1) {
                outfileOmega << t << "   " << k << "    " << Omega[0] << "    " << Omega[1] << "    " << Omega[2] << endl;
            }
        }
        outfileOmegaTot << t << "    " << a << "    " << OmegaTot[0] << "    " << OmegaTot[1] << "    " << OmegaTot[2] << endl;
    }
    outfileOmega.close();
    outfileOmegaTot.close();
    
    filename = "OmegaGWfin_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmegafin;
    outfileOmegafin.open(filename.c_str());
    
    // output the final GW spectrum
    a = at[Nt-1][1];
    H = Ht[Nt-1][1];
    if (expansion > 0) {
        Theta = 4.0*PI/(1.0*Nkhat)*3.0*pow(1.0/H,2.0)/(16.0*pow(PI*L,3.0));
    } else {
        Theta = 4.0*PI/(1.0*Nkhat)*3.0/(16.0*pow(PI*L,3.0));
    }
    for (int jk = 0; jk < Nk; jk++) {
        k = klist[jk];
        Omega = zeroNa;
        for (int ja = 0; ja < Na; ja++) {
            for (int jd = 0; jd < Nkhat; jd++) {
                for (int j6 = 0; j6 < 6; j6++) {
                    u0 = u[Nt-1][jd][jk][ja][j6];
                    du0 = du[Nt-1][jd][jk][ja][j6];
                    Omega[ja] += pow(k,3.0)*Theta*j6coef[j6]*(pow(abs(du0 + H*u0),2.0) + pow(abs(k*u0/a),2.0))/2.0;
                }
            }
        }
        outfileOmegafin << k << "    " << Omega[0] << "    " << Omega[1] << "    " << Omega[2] << endl;
    }
    outfileOmegafin.close();
    
    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    return 0;
}
