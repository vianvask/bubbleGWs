#include "functions.h"

int main (int argc, char *argv[]) {
    // nucleation rate parameters, units are chosen such that H0 = 1;
    const double beta = atof(argv[1]);
    const double gammapbeta = atof(argv[2]);
    const int index = atoi(argv[3]);
    
    cout << setprecision(5) << fixed;
    clock_t time_req = clock(); // timing
    rgen mt(time(NULL)*(index+1)/beta); // random number generator

    // bubble nucleation rate
    function<double(double)> Gamma = [beta, gammapbeta](double t) {
        return exp(beta*t - pow(gammapbeta*beta*t,2.0)/2.0);
    };
    
    // #nucleation sites
    int Nn = 100000;
    
    // #points on the bubble surfaces
    int Ns = 20000;
    
    // #timesteps
    const int Nt = 10000;
    
    // simulation volume determined by bar{N}(t=t_p) = J
    int J = 140;
    
    // largest GW wavenumber, Nk = k_max/k_min = #k values
    const int Nk = 240;
    
    // determine the time range and simulation volume
    vector<double> trange = findtrange(Gamma, 0.01, J, 6.0, 0.001);
    double tmin = trange[0];
    double tmax = trange[1];
    double dt = (tmax-tmin)/(1.0*Nt);
    double L = trange[2];
    double x1 = -L/2.0, x2 = L/2.0;
    
    cout << "time range: (" << tmin << ", " << tmax << ")" << endl;
    cout << "dt = " << dt << endl;
    cout << "L = " << L << endl;
    
    double tmaxnuc = min(trange[3],tmax);
    double dtnuc = (tmaxnuc-tmin)/(1.0*Nt);
    
    vector<vector<double> > Ft, taut, at, Ht, atau, ttau, Nb;
    
    // evolution of conformal time, scale factor, Hubble rate and expected number of bubbles, Nbar[jt][N,dN/dt]
    averageevolution(Gamma, tmin, Nt, dtnuc, Ft, taut, at, Ht, atau, ttau);
    Nb = Nbar(Gamma, x1, x2, Ft, taut, at);
    
    // nucleate bubbles
    vector<bubble> bubbles;
    int Ntry = 0;
    while ((bubbles.size() < 10 || bubbles.size() > J) && Ntry < 100) {
        bubbles = nucleate(x1, x2, Nn, taut, Nb, mt);
        cout << bubbles.size() << endl;
        Ntry++;
    }
    cout << "#bubbles = " << bubbles.size() << endl;
    
    // evolution of conformal time, scale factor and Hubble rate with larger dt
    Ft.clear(); taut.clear(); at.clear(); Ht.clear(); atau.clear(); ttau.clear();
    averageevolution(Gamma, tmin, Nt, dt, Ft, taut, at, Ht, atau, ttau);
    double taumax = taut[taut.size()-1][1];
    
    // generate a list of k values
    const double kmin = 1.0/(2.0*L);
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
    
    // initialize T_ij: Nt timesteps, Nkhat k directions, Nk k values, 3 approximations, 6 ij components
    int Na = 3;
    vector<vector<vector<vector<vector<complex<double> > > > > > T(Nt, vector<vector<vector<vector<complex<double> > > > > (Nkhat, vector<vector<vector<complex<double> > > > (Nk, vector<vector<complex<double> > > (Na, vector<complex<double> > (6, zero)))));
    
    // initialize h_ij and its time derivative in the same way as T_ij
    vector<vector<vector<vector<vector<complex<double> > > > > > u = T, du = T;
    
    string filename = "B_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(beta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileB;
    outfileB.open(filename.c_str());
    
    double tau, taun, tauc, t, tc, a, ac, H, R, Rc, kX, theta, phi;
    vector<double> xc(3, 0.0), xh(3, 0.0), X(3, 0.0), F(Na, 0.0), F6(6, 0.0);

    // loop over all the bubbles
    for (int jb = 0; jb < bubbles.size(); jb++) {
        cout << jb << endl;
        
        xc = bubbles[jb].x;
        taun = bubbles[jb].tau;
        
        // loop over the points on the sphere
        for (int jp = 0; jp < Ns; jp++) {
            xh = xhat[jp];
            
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
                outfileB << theta << "   " << phi << "   " << beta*Rc << "   " << ac << endl;
            }
            
            // compute the contribution to the stress-energy tensor
            for (int jt = 0; jt < Nt; jt++) {
                tau = taut[jt][1];
                t = taut[jt][0];
                if (tau > taun) {
                    a = at[jt][1];
                    R = radius(tau, taun);
                    for (int j=0; j<3; j++) {
                        X[j] = xc[j] + R*xh[j];
                    }
                    //X = periodic(x1, x2, X);
                                        
                    // 0: envelope, 1: xi = 2, 2: xi = 3
                    if (R < Rc) {
                        F[0] = 4.0*PI/(3.0*Ns)*pow(R,3.0);
                        F[1] = 4.0*PI/(3.0*Ns)*pow(R,3.0);
                        F[2] = 4.0*PI/(3.0*Ns)*pow(R,3.0);
                    } else {
                        F[0] = 0.0;
                        F[1] = 4.0*PI/(3.0*Ns)*pow(R,3.0)*pow(Rc/R,3.0);
                        F[2] = 4.0*PI/(3.0*Ns)*pow(R,3.0)*pow(Rc/R,4.0);
                    }
                    
                    for (int ja = 0; ja < Na; ja++) {
                        F6[0] = F[ja]*xh[0]*xh[0];
                        F6[1] = F[ja]*xh[0]*xh[1];
                        F6[2] = F[ja]*xh[0]*xh[2];
                        F6[3] = F[ja]*xh[1]*xh[1];
                        F6[4] = F[ja]*xh[1]*xh[2];
                        F6[5] = F[ja]*xh[2]*xh[2];
                        
                        for (int jk = 0; jk < Nk; jk++) {
                            k = klist[jk];
                            for (int jd = 0; jd < Nkhat; jd++) {
                                kX = k*inner(khat[jd],X);
                                for (int j6 = 0; j6 < 6; j6++) {
                                    T[jt][jd][jk][ja][j6] += F6[j6]*exp(-I*kX);
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
                        
                        u[jt+1][jd][jk][ja][j6] = (2.0*pow(dt,2.0)*(T0 - pow(k,2.0)*u0) + 2.0*pow(a,2.0)*(2.0*u0 - um1) + 3.0*pow(a,2.0)*dt*H*um1)/(pow(a,2.0)*(2.0 + 3.0*dt*H));
                    }
                }
            }
        }
    }
    
    // compute the TT projection of u_ij(k) and its time derivative
    for (int jt = 1; jt < Nt; jt++) {
        for (int jk = 0; jk < Nk; jk++) {
            for (int jd = 0; jd < Nkhat; jd++) {
                for (int ja = 0; ja < Na; ja++) {
                    u[jt][jd][jk][ja] = TTprojection6(u[jt][jd][jk][ja], khat[jd]);
                    
                    for (int j6 = 0; j6 < 6; j6++) {
                        du[jt][jd][jk][ja][j6] = (u[jt][jd][jk][ja][j6] - u[jt-1][jd][jk][ja][j6])/dt;
                    }
                }
            }
        }
    }
    
    filename = "OmegaGW_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(beta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmega;
    outfileOmega.open(filename.c_str());
    filename = "OmegaTotGW_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(beta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmegaTot;
    outfileOmegaTot.open(filename.c_str());
    
    const vector<double> zeroNa(Na, 0.0);
    const vector<double> j6coef {1.0, 2.0, 2.0, 1.0, 2.0, 1.0};
    vector<double> Omega(Na), OmegaTot(Na);
    double Theta;
    
    // output the GW spectrum in the end and the total GW energy density as a function of time
    for (int jt = 0; jt < Nt; jt++) {
        t = at[jt][0];
        a = at[jt][1];
        H = Ht[jt][1];
        Theta = 4.0*PI/(1.0*Nkhat)*3.0*pow(beta/H,2.0)/(16.0*pow(PI*L,3.0));
        
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
            
            if ((20*(jt+1))%Nt == 0) {
                outfileOmega << t << "   " << k/beta << "    " << Omega[0] << "    " << Omega[1] << "    " << Omega[2] << endl;
            }
        }
        outfileOmegaTot << t << "    " << a << "    " << OmegaTot[0] << "    " << OmegaTot[1] << "    " << OmegaTot[2] << endl;
    }
    outfileOmega.close();
    outfileOmegaTot.close();
    
    filename = "OmegaGWfin_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(beta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmegafin;
    outfileOmegafin.open(filename.c_str());
    
    a = at[Nt-1][1];
    H = Ht[Nt-1][1];
    Theta = 4.0*PI/(1.0*Nkhat)*3.0*pow(beta/H,2.0)/(16.0*pow(PI*L,3.0));
    for (int jk = 0; jk < Nk; jk++) {
        k = klist[jk];
        Omega = zeroNa;
        for (int ja = 0; ja < Na; ja++) {
            for (int jd = 0; jd < Nkhat; jd++) {
                for (int j6 = 0; j6 < 6; j6++) {
                    u0 = u[Nt-1][jd][jk][ja][j6];
                    du0 = du[Nt-1][jd][jk][ja][j6];
                    Omega[ja] += pow(k,3.0)*Theta*j6coef[j6]*0.5*(pow(abs(du0 + H*u0),2.0) + pow(abs(k*u0/a),2.0));
                }
            }
        }
        outfileOmegafin << k/beta << "    " << Omega[0] << "    " << Omega[1] << "    " << Omega[2] << endl;
    }
    outfileOmegafin.close();
    
    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    return 0;
}
