/*
 g++ -o hPT main.cpp basics.cpp functions.cpp -O3 -std=c++11
 ./hPT 8.0 0.0 0
*/

#include "functions.h"

int main (int argc, char *argv[]) {
    // nucleation rate parameters, units are chosen such that H0 = 1;
    const double beta = atof(argv[1]);
    const double gammapbeta = atof(argv[2]);
    const int index = atoi(argv[3]);
    
    clock_t time_req = clock();
    cout << setprecision(4) << fixed;

    rgen mt(time(NULL)*(1+beta)); // initialize random number generator
    
    // bubble nucleation rate
    function<double(double)> Gamma = [beta, gammapbeta](double t) {
        return exp(beta*t - pow(gammapbeta*beta*t,2.0)/2.0);
    };
    
    // number of bubbles including nucleation inside other bubbles
    int J = 50;
    
    // discretization of the bubble surfaces
    int Ns = 1000;
    
    // timesteps
    const int Nt = 3200;
    double tmin = -1.4; double tmax = 1.8;
    double dt = (tmax-tmin)/(1.0*Nt);
       
    // average evolution
    vector<vector<double> > Ft, taut, at, Ht, atau;
    vector<double> tmp(2);
    tmp = averageevolution(Gamma, tmin, Nt, dt, Ft, taut, at, Ht, atau);
    
    // simulation volume = cube of edge length L
    double L = 1.0/tmp[0];
    double x1 = -L/2.0, x2 = L/2.0;
    double lasttau = taut[taut.size()-1][1];
    
    // evolution of expected number of bubbles, Nbar[jt][N,dN/dt]
    vector<vector<double> > Nb;
    Nb = Nbar(Gamma, x1, x2, Ft, taut, at);
    
    bubble b0;
    vector<bubble> bubbles;
    double tn, taun, tauc, Rc, ac, a, tau, H, R, d, dj;
    vector<double> xh(3, 0.0), xc(3, 0.0), X(3, 0.0), F(2, 0.0), F6(6, 0.0);
    
    // try to nucleate J bubbles
    vector<int> jtj = jtlist(Nb, J, mt);
    bool flag;
    for (int j = 0; j < J; j++) {
        taun = taut[jtj[j]][1];
        for (int i = 0; i < 3; i++) {
            xc[i] = randomreal(x1,x2,mt);
        }
        flag = true;
        for (int jb = 0; jb < bubbles.size(); jb++) {
            d = distance(xc, bubbles[jb].x);
            R = radius(taun, bubbles[jb].tau);
            if (d < R) {
                flag = false;
            }
        }
        if (flag) {
            b0.x = xc;
            b0.tau = taun;
            bubbles.push_back(b0);
        }
    }
    J = bubbles.size();
    cout << "#bubbles = " << J << endl;
    
    // generate a list of k values
    const double kmin = 1.0/(2.0*L);
    double k = kmin;
    vector<double> klist;
    while (k < 100.1*kmin) {
        klist.push_back(k);
        k += kmin;
    }
    const int Nk = klist.size();
    
    // generate points on a unit sphere
    vector<vector<double> > xhcoll = sphereN(Ns);
    Ns = xhcoll.size();
    
    // initialize T_ij: Nt timesteps, 3 k directions, Nk k values, 2 approximations, 6 independent components
    vector<vector<vector<vector<vector<complex<double> > > > > > T(Nt, vector<vector<vector<vector<complex<double> > > > > (3, vector<vector<vector<complex<double> > > > (Nk, vector<vector<complex<double> > > (2, vector<complex<double> > (6, zero)))));
    
    // initialize h_ij and its time derivative in the same way as T_ij
    vector<vector<vector<vector<vector<complex<double> > > > > > u, du;
    u = T; du = T;

    // loop over all the bubbles
    for (int jb = 0; jb < bubbles.size(); jb++) {
        cout << jb << endl;
        
        xc = bubbles[jb].x;
        taun = bubbles[jb].tau;
        
        // loop over the points on the sphere
        for (int jp = 0; jp < Ns; jp++) {
            xh = xhcoll[jp];
            
            // find when the collision happens in the direction xh
            tauc = findtauc(x1, x2, taun, xh, xc, bubbles, jb, 2.0*lasttau);
            ac = interpolate(tauc, atau);
            Rc = radius(tauc, taun);
            
            // compute the contribution to the stress-energy tensor
            for (int jt = 0; jt < Nt; jt++) {
                tau = taut[jt][1];
                if (tau > taun) {
                    a = at[jt][1];
                    R = radius(tau, taun);
                    for (int j=0; j<3; j++) {
                        X[j] = xc[j] + R*xh[j];
                    }
                    //X = periodic(x1, x2, X);
                                        
                    // 0: envelope, 1: xi = 3
                    if (R < Rc) {
                        F[0] = 4.0*PI/(3.0*Ns)*pow(a*R,3.0);
                        F[1] = 4.0*PI/(3.0*Ns)*pow(a*R,3.0);
                    } else {
                        F[0] = 0.0;
                        F[1] = 4.0*PI/(3.0*Ns)*pow(ac*Rc,3.0)*(ac*Rc)/(a*R);
                    }
                    
                    for (int ja=0; ja<2; ja++) {
                        F6[0] = F[ja]*xh[0]*xh[0];
                        F6[1] = F[ja]*xh[0]*xh[1];
                        F6[2] = F[ja]*xh[0]*xh[2];
                        F6[3] = F[ja]*xh[1]*xh[1];
                        F6[4] = F[ja]*xh[1]*xh[2];
                        F6[5] = F[ja]*xh[2]*xh[2];
                        for (int jk = 0; jk < Nk; jk++) {
                            k = klist[jk];
                            for (int jd = 0; jd < 3; jd++) {
                                for (int j6 = 0; j6 < 6; j6++) {
                                    T[jt][jd][jk][ja][j6] += F6[j6]*exp(-I*k*X[jd]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // compute u_ij(k) and its time derivative
    double a1, a2, H1, H2;
    complex<double> k1, k2, l1, l2;
    for (int jt = 0; jt < Nt-1; jt++) {
        a1 = at[jt][1];
        H1 = Ht[jt][1];
        a2 = at[jt+1][1];
        H2 = Ht[jt+1][1];
        
        for (int jk = 0; jk < Nk; jk++) {
            k = klist[jk];
            for (int jd = 0; jd < 3; jd++) {
                for (int ja = 0; ja < 2; ja++) {
                    for (int j6 = 0; j6 < 6; j6++) {
                        k1 = du[jt][jd][jk][ja][j6];
                        l1 = T[jt][jd][jk][ja][j6]/(a1*a1) - 3.0*H1*k1 - pow(k/a1,2.0)*u[jt][jd][jk][ja][j6];
                        
                        k2 = du[jt][jd][jk][ja][j6] + dt*l1;
                        u[jt+1][jd][jk][ja][j6] = u[jt][jd][jk][ja][j6] + dt*(k1+k2)/2.0;
                        
                        l2 = T[jt+1][jd][jk][ja][j6]/(a2*a2) - 3.0*H2*k2 - pow(k/a2,2.0)*u[jt+1][jd][jk][ja][j6];
                        du[jt+1][jd][jk][ja][j6] = du[jt][jd][jk][ja][j6] + dt*(l1+l2)/2.0;
                    }
                }
            }
        }
    }
    
    // unit wavevectors
    vector<vector<double> > khat(3, vector<double> (3, 0.0));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            khat[i][j] = delta(i,j);
        }
    }
        
    // compute the TT projection of u_ij(k) and its time derivative
    for (int jt = 0; jt < Nt; jt++) {
        for (int jk = 0; jk < Nk; jk++) {
            for (int jd = 0; jd < 3; jd++) {
                for (int ja = 0; ja < 2; ja++) {
                    u[jt][jd][jk][ja] = TTprojection6(u[jt][jd][jk][ja], khat[jd]);
                    du[jt][jd][jk][ja] = TTprojection6(du[jt][jd][jk][ja], khat[jd]);
                }
            }
        }
    }
    
    // output the GW spectrum in the end and the total GW energy density as a function of time
    string filename = "OmegaGWbeta" + to_string_prec(beta,2) + "j" + to_string(index) + ".dat";
    ofstream outfileOmega;
    outfileOmega.open(filename.c_str());
    filename = "OmegaTotGWbeta" + to_string_prec(beta,2) + "j" + to_string(index) + ".dat";
    ofstream outfileOmegaTot;
    outfileOmegaTot.open(filename.c_str());
    
    vector<double> Omega(2), OmegaTot(2), k2Pu(2), Pdu(2), zero2(2,0.0);
    double Theta;
    for (int jt = 0; jt < Nt; jt++) {
        a = at[jt][1];
        H = at[jt][1];
        Theta = 4.0*PI/3.0*pow(H,-2.0)*3.0*pow(beta,2.0)/(16.0*pow(PI*L,3.0));
        
        OmegaTot = zero2;
        for (int jk = 0; jk < Nk; jk++) {
            k = klist[jk];
            Pdu = zero2;
            k2Pu = zero2;
            for (int ja = 0; ja < 2; ja++) {
                for (int jd = 0; jd < 3; jd++) {
                    for (int j6 = 0; j6 < 6; j6++) {
                        Pdu[ja] += (1.0 + ((j6==1) || (j6==2) || (j6==4)))*pow(abs(du[jt][jd][jk][ja][j6]),2.0);
                        k2Pu[ja] += (1.0 + ((j6==1) || (j6==2) || (j6==4)))*pow(k/a*abs(u[jt][jd][jk][ja][j6]),2.0);
                    }
                }
                OmegaTot[ja] += pow(k,2.0)*kmin*Theta*(Pdu[ja] + k2Pu[ja]);
            }
            if (jt == Nt - 1) {
                for (int ja = 0; ja < 2; ja++) {
                    Omega[ja] = pow(k,3.0)*Theta*(Pdu[ja] + k2Pu[ja]);
                }
                outfileOmega << k/beta << "    " << Omega[0] << "    " << Omega[1] << endl;
            }
        }
        outfileOmegaTot << at[jt][0] << "    " << a << "    " << OmegaTot[0] << "    " << OmegaTot[1] << endl;
    }
    outfileOmega.close();
    outfileOmegaTot.close();
    
    time_req = clock() - time_req;
    cout << "total evaluation time: " << ((double) time_req/CLOCKS_PER_SEC/60.0) << " minutes." << endl;
    return 0;
}
