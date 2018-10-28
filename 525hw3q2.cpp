 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include "E:\study\newmat10\newmatap.h"           
#include "E:\study\newmat10\newmat.h"
#include "E:\study\newmat10\newmatio.h"
using namespace std;

double xmin = -2.0; double xmax = 3; double t = 0; double tmax = 1; double dt = 0.01;

//solver(0.05);
//solver(0.01);
//solver(0.005);

float max(float a, float b) {
	return (b < a) ? a : b;
}

ColumnVector maxmatrix(ColumnVector vect, double var1) {
	for (int i = 1; i <= vect.Nrows(); i++) {
		if (vect(i) < var1) {
			vect(i) = var1;
		}
	}
	return vect;
}


double function(double inputx) {
	if (inputx < -1) { return 0.0;}
	if (inputx <= 0) { return inputx + 1.0; }
	return  1.0; 
}

Matrix solver (double dx) {
	int N = (int) (5.0 / dx);
	int M = 100;
	double k1 = 1 - dt / dx; double k2 = dt / dx;
	
	Matrix solution (N+1,M+1);// This is a N+1 x M+1 matrix
	for (int i = 1; i <= N + 1; i++) {
		for (int j = 1; j <= M + 1; j++) {
			solution(i, j) = 0;
		}
	}
	vector <double> vetx; vetx.push_back(0);
	for (int i = 0; i <= N+1; i++) {
		vetx.push_back( xmin + i*dx);
		//cout << vetx[i] << "  ";
	}
	for (int i = 1; i <= N + 1; i++) {
		solution(i,1) = function(vetx[i]);
	}
	double fixed = solution(1,1);
	//for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= M; j++) {
			ColumnVector temp (1);
			temp(1) = fixed;
			Matrix sub = temp & (solution.Column(j)).Rows(1,N);
			solution.Column(j + 1) = k1*solution.Column(j )+ k2*sub;

		}
	//}

	return solution;
}

double baroption(double s0, double K, double r, double T, double sigma,
	double sb, double ds) {
	double dt1 = 1.0 / 1200.0; ds = 0.5;
	double smax = 100.0;
	int M = round((smax - sb) / ds);
	ds = (smax - sb) / M;
	int N = round(T / dt1);
	dt1 = T / N;
	Matrix matval(M + 1, N + 1);  matval = 0.0;
	ColumnVector vets(M + 1);
	for (int i = 1; i <= M + 1; i++) {
		//vets(i) = sb + (i-1)* (smax - sb) / (M + 1);
		vets(i) = sb + (i - 1)* (0.5);
	}
	ColumnVector veti = vets / ds;
	ColumnVector vetj(N + 1);
	for (int i = 1; i <= N + 1; i++) {
		vetj(i) = i - 1;

	}
	// set up boundary
	matval.Column(N + 1) = maxmatrix(K - vets, 0);
	matval.Row(1) = 0;
	matval.Row(M + 1) = 0;
	//set up coef Matrix
	ColumnVector alpha = 0.25*dt1*(sigma*sigma*SP(veti, veti) - r*veti);
	ColumnVector beta = -dt1*0.5*(sigma*sigma*SP(veti, veti) + r);
	ColumnVector gamma = 0.25*dt1*(sigma*sigma*SP(veti, veti) + r*veti);
	//beta = beta.Rows(2, M);
	//alpha = alpha.Rows(2, M);
	//gamma = gamma.Rows(2, M);
	
	Matrix M1(M - 1, M - 1); M1 = 0.0;
	for (int i = 1; i <= M - 1; i++) {
		M1(i, i) = 1 - beta(i+1);
		if (i > 1) {
			M1(i, i - 1) = -alpha(i+1);
			M1(i - 1, i) = -gamma(i);
		}
		//for (int j = 1; j <= M - 1; j++)
			//M1(i, i) = 1 - beta(i);
	}
	// -(alpha.Rows(3,M), -1)+ (1 - beta.Rows(2, M)) - (gamma.Rows(2,M - 1), 1);
	Matrix M2(M - 1, M - 1); M2 = 0.0;
	for (int i = 1; i <= M - 1; i++) {
		M2(i, i) = 1 + beta(i+1);
		if (i > 1) {
			M2(i, i - 1) = alpha(i+1);
			M2(i - 1, i) = gamma(i);
		}
		//diag(alpha(3:M), -1) + diag(l + beta(2:M)) + diag(gamma(2:M - 1), 1);
	}
	for (int j = N; j >= 1; j--) {
		//matval.Rows(2, M).Column(j)
		ColumnVector temp = (M1.i()*(M2*(matval.Rows(2, M).Column(j + 1))));
		//cout << temp.Nrows() << endl;
		for (int i = 1; i <= temp.Nrows(); i++) {
			matval(i+1,j)= temp(i);

		}
		//matval.Rows(2, M).Column(j);
		//cout << temp2.Nrows() << endl;
	}
	//Matrix temp = M1.i()*M2;
	//matval.Column(500).Rows(2, M) = M1.i()*(M2*(matval.Column(500 + 1).Rows(2, M)));
	//double theprice = matval(21,1);
	//cout << matval.Column(500 + 1).Rows(2, M) << endl;
	//cout << M2*matval.Column(500 + 1).Rows(2, M) << endl;
	//cout << M2(119,119) << endl; cout << M2(118, 119) << endl; cout << M2(119, 118) << endl;
	//cout << (matval.Column(20).Rows(2, M)).Nrows() << endl;
	//cout << vets << endl;
	return matval(21,1);
}

double normcdf(const double& z) {
	if (z > 8.0) { return 1.0; }; // this guards against overflow 
	if (z < -8.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double formula(double S0, double K, double r, double T, double sigma,
	double Sb) {
	double s2 = pow(sigma, 2);
	double a = pow((Sb / S0) , (-1 + (2 * r / s2)));
	double b = pow((Sb / S0) , (1 + (2 * r / s2)));
	double d1 = (log(S0 / K) + (r + s2 / 2)* T) / (sigma*sqrt(T)) ;
	double d2 = (log(S0 / K) + (r - s2 / 2)* T) / (sigma*sqrt(T));
	double d3 = (log(S0 / Sb) + (r + s2 / 2)* T) / (sigma*sqrt(T));
	double d4 = (log(S0 / Sb) + (r - s2 / 2)* T) / (sigma*sqrt(T));
	double d5 = (log(S0 / Sb) - (r - s2 / 2)* T) / (sigma*sqrt(T));
	double d6 = (log(S0 / Sb) - (r + s2 / 2)* T) / (sigma*sqrt(T));
	double d7 = (log(S0*K / (Sb*Sb)) - (r - s2 / 2)* T) / (sigma*sqrt(T));
	double d8 = (log(S0*K / (Sb*Sb)) - (r + s2 / 2)* T) / (sigma*sqrt(T));
	double P = K*exp(-r*T) * (normcdf(d4) - normcdf(d2) - a*(normcdf(d7) - normcdf(d5)))
		- S0* (normcdf(d3) - normcdf(d1) - b* (normcdf(d8) - normcdf(d6)));
	return P;


}

int main(int argc, char* argv[]) {
	
	//uncomment for Q2
	/*Matrix sol1 = solver(0.05);
	Matrix sol2 = solver(0.01);
	Matrix sol3 = solver(0.005);
	ofstream outfile1("dx005");
	ofstream outfile2("dx001");
	ofstream outfile3("dx0005");
	outfile1 <<  setprecision(10) << sol1;
	outfile2 <<  setprecision(10) << sol2;
	outfile3 <<  setprecision(10) << sol3;*/

	//Q3
	double barprice = baroption(50,50,0.1,5.0/12.0,0.4,40,0.5);
	cout << barprice << endl;
	double formulaprice = formula(50.0, 50.0, 0.1, 5.0 / 12.0, 0.4, 40);
	cout << formulaprice << endl;



}


