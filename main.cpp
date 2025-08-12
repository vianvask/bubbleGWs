#include "functions.h"

int main (int argc, char *argv[]) {
    
    // nucleation rate parameters
    const double beta = atof(argv[1]);
    const double gammapbeta = atof(argv[2]);
    // file naming
    const int index = atoi(argv[3]);
    // 0: static Minkowski space, 1: FLRW space with transition from vacuum to radiation dominance
    const int expansion = atoi(argv[4]);
    
    cout << setprecision(5) << fixed;
    clock_t time_req = clock(); // timing
    rgen mt(time(NULL)+(index+1)+beta); // random number generator

    string filename;
    ofstream outfile;
    
    // bubble nucleation rate, units chosen such that H0 = 1
    function<double(double)> Gamma = [beta, gammapbeta](double t) {
        return exp(beta*t - pow(gammapbeta*beta*t,2.0)/2.0);
    };
    
    int Nn = 10000; // #nucleation sites
    int Ns = 10000; // #points on the bubble surfaces
    int Nt = 5000; // #timesteps
    int Nk = 50; // #k values
    
    int J = 20; // bar{N}(t=t_p) = J, fixes L
    double barNtmin = 0.001; // bar{N}(t=t_min) = barNtmin, fixes t_min
    double ftmax = 50.0; // t_max = t_p + ftmax*<R>
    double barFtmaxnuc = 0.01; // bar{F}(t=t_max,nuc) = barFtmaxnuc, fixes t_max,nuc

    // determine the time range and simulation volume
    vector<double> trange = findtrange(Gamma, barNtmin, J, ftmax, barFtmaxnuc, expansion);
    double tmin = trange[0];
    double tmax = trange[1];
    double dt = (tmax-tmin)/(1.0*Nt);
    double tmaxnuc = min(trange[3],tmax);
    double dtnuc = (tmaxnuc-tmin)/(1.0*Nt);
    double L = trange[2];
    double x1 = -L/2.0, x2 = L/2.0;
        
    // evolution of conformal time, scale factor and Hubble rate
    vector<vector<double> > Ft, taut, at, Ht, ttau;
    averageevolution(Gamma, tmin, Nt, dtnuc, Ft, taut, at, Ht, ttau, expansion);
    
    if (index == 0) {
        filename = "P_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + ".dat";
        writeToFile(Ft, filename);
        
        filename = "a_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + ".dat";
        writeToFile(at, filename);
        
        filename = "tau_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + ".dat";
        writeToFile(taut, filename);
        
        filename = "H_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + ".dat";
        writeToFile(Ht, filename);
        
        vector<vector<double> > pRc = RcPDF(Gamma, Ft, taut, at);
        filename = "pRc_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + ".dat";
        writeToFile(pRc, filename);
    }
        
    // nucleate bubbles
    vector<bubble> bubbles;
    int Ntry = 0;
    while ((bubbles.size() < J/2.0 || bubbles.size() > 2.0*J) && Ntry < 100) {
        bubbles = nucleate(Gamma, x1, x2, Nn, taut, at, mt);
        cout << bubbles.size() << endl;
        Ntry++;
    }
    cout << "#bubbles = " << bubbles.size() << endl;
    
    // recompute the background evolution with up to larger times with a larger dt
    Ft.clear(); taut.clear(); at.clear(); Ht.clear(); ttau.clear();
    averageevolution(Gamma, tmin, Nt, dt, Ft, taut, at, Ht, ttau, expansion);
    double taumax = taut[taut.size()-1][1];
    
    // generate a list of k values in log scale
    double kmin = 1.0/L;
    double kmax = 30.0*beta;
    double dlogk = (log(kmax) - log(kmin))/(1.0*(Nk-1));
    double k = kmin;
    vector<double> klist;
    for (int jk = 0; jk < Nk; jk++) {
        klist.push_back(k);
        k = exp(log(k) + dlogk);
    }
    
    // generate x and k directions
    vector<vector<double> > xhat = sphereN(Ns);
    int Nkhat = 3;
    vector<vector<double> > khat(Nkhat, vector<double> (3, 0.0));
    for (int jd = 0; jd < Nkhat; jd++) {
        for (int j = 0; j < 3; j++) {
            khat[jd][j] = delta(jd,j);
        }
    }
    
    // initialize T_ij: Nt timesteps, 6 k directions, Nk k values, Na approximations, 6 ij components
    int Na = 3;
    vector<vector<vector<vector<vector<complex<double> > > > > > T(Nt, vector<vector<vector<vector<complex<double> > > > > (Nkhat, vector<vector<vector<complex<double> > > > (Nk, vector<vector<complex<double> > > (Na, vector<complex<double> > (6, zero)))));
    
    // file for collision times on the surface of the first bubble
    filename = "B_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    outfile.open(filename.c_str());
    
    // initialize binning of the collision radius
    const int Nbins = 1000;
    const double Rcmax = 20.0/beta, acRcmax = 20.0/beta;
    vector<double> Rcbins(Nbins, 0.0), acRcbins(Nbins, 0.0);
    
    double tau, taun, tauc, t, tc, a, an, ac, H, Hc, R, dR, Rc, dA, theta, phi, K;
    complex<double> eikX;
    vector<double> xc(3, 0.0), xh(3, 0.0), X(3, 0.0), xixj(6, 0.0), F(Na, 0.0);

    // loop over all the bubbles
    for (int jb = 0; jb < bubbles.size(); jb++) {
        cout << jb << endl;
        
        xc = bubbles[jb].x;
        taun = bubbles[jb].tau;
        an = bubbles[jb].a;
        
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
            ac = interpolate(tc, at);
            Hc = interpolate(tc, Ht);
            Rc = radius(tauc, taun);
            
            // output the collision radii of first bubble
            if (jb < 1) {
                theta = acos(xh[2]);
                phi = sgn(xh[1]);
                if (sqrt(xh[0]*xh[0]+xh[1]*xh[1]) > 0) {
                    phi = sgn(xh[1])*acos(xh[0]/sqrt(xh[0]*xh[0]+xh[1]*xh[1]));
                }
                outfile << theta << "   " << phi << "   " << Rc << "   " << ac << endl;
            }
            
            // add the point to the distribution of collision radii
            if (Rc < Rcmax) {
                Rcbins[(int) floor(Nbins*Rc/Rcmax)] += pow(ac,4.0);
            }
            if (ac*Rc < acRcmax) {
                acRcbins[(int) floor(Nbins*ac*Rc/acRcmax)] += pow(ac,4.0);
            }
            
            // compute the contribution to the stress-energy tensor
            R = 0.0;
            K = 0.0;
            for (int jt = 0; jt < Nt; jt++) {
                tau = taut[jt][1];
                if (tau > taun) {
                    a = at[jt][1];
                    H = Ht[jt][1];
                    dR = radius(tau, taun) - R;
                    R = radius(tau, taun);
                    for (int j = 0; j < 3; j++) {
                        X[j] = xc[j] + R*xh[j];
                    }
                    
                    // dA = dOmega (a R)^2/a
                    dA = 4.0*PI/(a*Ns)*pow(a*R,2.0);
                                        
                    // K = gamma*sigma/DeltaV, solve dK/dR + 2 K/R (1 + 3/2 a R H) = a or 0
                    if (R < Rc) {
                        if (expansion == 0) {
                            K = a*R/3.0;
                        } else {
                            K += (a - 2.0*K/R*(1.0+3.0/2.0*a*H*R))*dR;
                        }
                        F[0] = K;
                        F[1] = K;
                        F[2] = K;
                    } else {
                        if (expansion == 0) {
                            K = pow(Rc/R,2.0)*ac*Rc/3.0;
                        } else {
                            K += -2.0*K/R*(1.0+3.0/2.0*a*H*R)*dR;
                        }
                        F[0] = 0.0; // envelope approximation
                        F[1] = K; // bulk flow approximation
                        F[2] = K*ac*Rc/(a*R); // extra a_c*R_c/a*R dissipation
                    }
                    
                    for (int jk = 0; jk < Nk; jk++) {
                        k = klist[jk];
                        for (int jd = 0; jd < Nkhat; jd++) {
                            eikX = exp(-I*k*inner(khat[jd],X));
                            for (int ja = 0; ja < Na; ja++) {
                                for (int j6 = 0; j6 < 6; j6++) {
                                    T[jt][jd][jk][ja][j6] += dA*F[ja]*xixj[j6]*eikX;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    outfile.close();
    
    // output the R_c and a_c R_c distributions
    ofstream outfileRc, outfileacRc;
    filename = "Rc_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    outfileRc.open(filename.c_str());
    filename = "acRc_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    outfileacRc.open(filename.c_str());
    for (int jb = 0; jb < Nbins; jb++) {
        outfileRc << beta*jb*Rcmax/(1.0*Nbins) << "    " << Rcbins[jb] << endl;
        outfileacRc << beta*jb*acRcmax/(1.0*Nbins) << "    " << acRcbins[jb] << endl;
    }
    outfileRc.close();
    outfileacRc.close();
    
    // TT projection
    for (int jt = 0; jt < Nt; jt++) {
        for (int jk = 0; jk < Nk; jk++) {
            for (int jd = 0; jd < Nkhat; jd++) {
                for (int ja = 0; ja < Na; ja++) {
                    T[jt][jd][jk][ja] = TTprojection6(T[jt][jd][jk][ja], khat[jd]);
                }
            }
        }
    }
    
    // initialize u and du, they have half the timesteps of T because RK4 method is used
    int jtu = 0;
    for (int jt = 0; jt < Nt-2; jt+=2) {
        jtu++;
    }
    vector<vector<vector<vector<vector<complex<double> > > > > > u(jtu, vector<vector<vector<vector<complex<double> > > > > (Nkhat, vector<vector<vector<complex<double> > > > (Nk, vector<vector<complex<double> > > (Na, vector<complex<double> > (6, zero)))));
    vector<vector<vector<vector<vector<complex<double> > > > > > du = u;
    
    // solve the perturbed Einstein equation
    vector<vector<complex<double> > > Tt(Nt, vector<complex<double> > (2,zero));
    for (int jt = 0; jt < Nt; jt++) {
        Tt[jt][0] = at[jt][0];
    }
    vector<complex<double> > udu(2, zero);
    for (int jk = 0; jk < Nk; jk++) {
        k = klist[jk];
        for (int jd = 0; jd < Nkhat; jd++) {
            for (int ja = 0; ja < Na; ja++) {
                for (int j6 = 0; j6 < 6; j6++) {
                    for (int jt = 0; jt < Nt; jt++) {
                        Tt[jt][1] = T[jt][jd][jk][ja][j6];
                    }
                    jtu = 0;
                    udu[0] = zero;
                    udu[1] = zero;
                    for (int jt = 0; jt < Nt-2; jt+=2) {
                        udu = PEEstep(jt, udu, Tt, at, Ht, k);
                        u[jtu][jd][jk][ja][j6] = udu[0];
                        du[jtu][jd][jk][ja][j6] = udu[1];
                        jtu++;
                    }
                }
            }
        }
    }
    Nt = jtu; // Nt to match the length of u and du vectors
    
    // file for the GW spectrum as a function of time
    filename = "OmegaGW_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmega;
    outfileOmega.open(filename.c_str());
    
    // file for the total GW abundance as a function of time
    filename = "OmegaTotGW_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmegaTot;
    outfileOmegaTot.open(filename.c_str());
    
    const vector<double> zeroNa(Na, 0.0);
    const vector<double> j6coef {1.0, 2.0, 2.0, 1.0, 2.0, 1.0};
    complex<double> u0, du0;
    double Theta, dtheta = 4.0*PI/Nkhat;
    vector<double> Omega(Na), OmegaTot(Na);
    
    // output the GW spectrum and the total GW energy density as a function of time
    for (int jt = 0; jt < Nt; jt++) {
        t = at[2*jt][0];
        a = at[2*jt][1];
        H = Ht[2*jt][1];
        if (expansion > 0) {
            Theta = dtheta*3.0*pow(1.0/H,2.0)/(16.0*pow(PI*L,3.0));
        } else {
            Theta = dtheta*3.0/(16.0*pow(PI*L,3.0));
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
                    OmegaTot[ja] += dlogk*Omega[ja];
                }
            }
            outfileOmega << t << "   " << k/beta << "    " << Omega[0] << "    " << Omega[1] << "    " << Omega[2] << endl;
        }
        outfileOmegaTot << t << "    " << a << "    " << OmegaTot[0] << "    " << OmegaTot[1] << "    " << OmegaTot[2] << endl;
    }
    outfileOmega.close();
    outfileOmegaTot.close();
    
    // file for the final GW spectrum
    filename = "OmegaGWfin_beta_" + to_string_prec(beta,2) + "_gammaperbeta_" + to_string_prec(gammapbeta,2) + "_j_" + to_string(index) + ".dat";
    ofstream outfileOmegafin;
    outfileOmegafin.open(filename.c_str());
    
    // output the final GW spectrum
    a = at[2*(Nt-1)][1];
    H = Ht[2*(Nt-1)][1];
    if (expansion > 0) {
        Theta = dtheta*3.0*pow(1.0/H,2.0)/(16.0*pow(PI*L,3.0));
    } else {
        Theta = dtheta*3.0/(16.0*pow(PI*L,3.0));
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
        outfileOmegafin << k/beta << "    " << Omega[0] << "    " << Omega[1] << "    " << Omega[2] << endl;
    }
    outfileOmegafin.close();
     
    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    return 0;
}
