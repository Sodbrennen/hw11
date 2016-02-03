#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void CrankStep(cmplx* const psi1, cmplx* const psi0, double* V, const double m, const double hq, const double dt, const double dx, const int Nx, const double xmin);
//-----------------------------------

int main() {

	const int Nx = 300;
	const double xmin = -40;
	const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax - xmin)/(Nx - 1);
	const double dt = dx/100.0;
	double t = 0;
 	const int Na = 10;
 	int Nk = int(Tend / Na / dt + 0.5);
 
 	const double lambda = 10;
	const double omega = 0.2;
	const double m = 1;
	const double hq = 1;
	const double k = omega*omega*m;
	const double alpha = pow(m*k/hq, 1/4.0);

	double* V = new double[Nx];
	
	for(int i = 0; i < Nx; i++){
		V[i] = 0.5*k*(xmin + i*dx)*(xmin + i*dx);
 	}

	stringstream strm;
 
 	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);
 
	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);
	writeToFile(psi0, "psi_0", dx, Nx, xmin, alpha, lambda, omega, t);
 
 
 
 	for (int i = 1; i <= Na; i++) {
 		for (int j = 1; j <= Nk-1; j++) {
		  
		    CrankStep(psi1, psi0, V, m, hq, dt, dx, Nx, xmin);
		    
		    h = psi0;
		    psi0 = psi1;
		    psi1 = h;
 
         t+=dt;
		    t += dt;
 		}
		
 		strm.str("");
 		strm << "psi_" << i;
 		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
 	}
  cout << "t = " << t << endl;
	
	cout << "t = " << t << endl;
	
	delete[] psi0;
	delete[] psi1;
	delete[] V;
 
 	return 0;
}
//-----------------------------------
void CrankStep(cmplx* const psi1, cmplx* const psi0, double* V, const double m, const double hq, const double dt, const double dx, const int Nx, const double xmin) {

	double x;
	cmplx* d = new cmplx[Nx];
	cmplx* r = new cmplx[Nx];
	cmplx a = - cmplx(0,1)*hq*dt/(4*m*dx*dx);

	for(int i = 0; i < Nx; i++) 
		d[i] = 1.0 + cmplx(0,1)*(hq*dt/(2*m*dx*dx) + dt*V[i]/(2*hq));
	

	r[0] = psi0[0] - a*(psi0[1] - cmplx(2,0)*psi0[0]) - cmplx(0,1)*dt/2.0*V[0]*psi0[0];	

	for(int i = 1; i < Nx - 1; i++) {
 	  	x = xmin + i * dx;
	  	r[i] = psi0[i] - a*(psi0[i+1] - cmplx(2,0)*psi0[i] + psi0[i-1]) - cmplx(0,1)*dt/2.0*V[i]*psi0[i];
	}
	
	r[Nx-1] = psi0[Nx-1] - a*(-cmplx(2,0)*psi0[Nx-1] + psi0[Nx-2]) - cmplx(0,1)*dt/2.0*V[Nx-1]*psi0[Nx-1];


	
	for(int i = 1; i < Nx; i++){
		d[i] = d[i] - a*a/d[i-1];		
		r[i] = r[i] - a*r[i-1]/d[i-1];
	}
	
	psi1[Nx-1] = r[Nx-1]/d[Nx-1];
	
	for(int i = Nx-2; i > 0; i--){
		psi1[i] = (r[i] - a*psi1[i+1])/d[i];
	}
	     

	delete[] d;
	delete[] r;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
