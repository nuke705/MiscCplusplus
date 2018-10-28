#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
using namespace std;

float u, p, p_u, p_d, r, K;
float d, downtick_prob;
float S10,S20, T, volatility, target = 9.472992;
int N;
float max(float a, float b) {
	return (b < a) ? a : b;
}

double euro_spread_price(int N) {
	double  S10 = 101, S20 = 100, K = 3.0, volatility1 = 0.4, volatility2 = 0.2;
	double r = 0.05, T = 0.5, deltaT = T / N, corr = 0.1;

	double u1 = exp(volatility1*sqrt(deltaT)), d1 = 1 / u1, p1 = (exp(r*deltaT) - d1) / (u1 - d1);
	double discount = exp(-r*deltaT), p_u1 = discount*p1, p_d1 = discount*(1 - p1);

	double u2 = exp(volatility2*sqrt(deltaT)), d2 = 1 / u2, p2 = (exp(r*deltaT) - d2) / (u2 - d2);

	double v1 = r - volatility1*volatility1 / 2, v2 = r - volatility2*volatility2 / 2;
	double puu = (0.25)*(1 + sqrt(deltaT)*(v1 / volatility1 + v2 / volatility2) + corr);
	double pud = (0.25)*(1 + sqrt(deltaT)*(v1 / volatility1 - v2 / volatility2) - corr);
	double pdu = (0.25)*(1 + sqrt(deltaT)*(-v1 / volatility1 + v2 / volatility2) - corr);
	double pdd = (0.25)*(1 + sqrt(deltaT)*(-v1 / volatility1 - v2 / volatility2) + corr);

	vector <double> S1(2 * N + 2); vector <double> S2(2 * N + 2);
	S1[1] = S10*pow(d1, N); S2[1] = S20*pow(d2, N);
	for (int i = 2; i <= 2 * N + 1; i++) {
		S1[i] = u1*S1[i - 1]; S2[i] = u2*S2[i - 1];
	}
	
	vector <vector <double>> call;// This is a N+1 x N+1 matrix
	for (int i = 0; i <= N * 2 + 1; i++) {
		vector <double> temp;
		call.push_back(temp);
		for (int j = 0; j <= N * 2 + 1; j++) {
			call[i].push_back(0);
		}
	}
	for (int i = 1; i <= 2 * N + 1; i += 2) {
		for (int j = 1; j <= 2 * N + 1; j += 2) {
			call[i][j] = max(S1[i] - S2[j] - K, 0);	
		}
	}
	
	for (int tau = 1; tau <= N; tau++) {
		for (int i = (tau + 1); i <= (2 * N + 1 - tau); i += 2) {
			for (int j = (tau + 1); j <= (2 * N + 1 - tau); j += 2) {
				double hold = puu*(call[i + 1][j + 1]) + pud*(call[i + 1][j - 1]) + 
					pdu*(call[i - 1][j + 1]) + pdd*(call[i - 1][j - 1]);
				call[i][j] = max( S1[i] - S2[j] - K, hold);
			}
		}
	}
	double price = call[N+1][N+1];
	return price;

}
int main(int argc, char* argv[]) {
	ofstream outfile("spreadcall");
	ofstream outfile2("timespent");
	//double price15 = euro_spread_price(15);
	clock_t startTime1 = clock();
	outfile << euro_spread_price(15) << endl;
	clock_t endTime1 = clock();
	clock_t clockTicksTaken1 = endTime1 - startTime1;
	outfile2<<  clockTicksTaken1 / (double)CLOCKS_PER_SEC << endl;

	clock_t startTime2 = clock();
	outfile << euro_spread_price(30) << endl;
	clock_t endTime2 = clock();
	clock_t clockTicksTaken2 = endTime2 - startTime2;
	outfile2 << clockTicksTaken2 / (double)CLOCKS_PER_SEC << endl;

	clock_t startTime = clock();
	outfile << euro_spread_price(60) << endl;
	clock_t endTime = clock();
	clock_t clockTicksTaken = endTime - startTime;
	outfile2 << clockTicksTaken / (double)CLOCKS_PER_SEC << endl;

	clock_t startTime4 = clock();
	outfile << euro_spread_price(120) << endl;
	clock_t endTime4 = clock();
	clock_t clockTicksTaken4 = endTime4 - startTime4;
	outfile2 << clockTicksTaken4 / (double)CLOCKS_PER_SEC << endl;

	clock_t startTime5 = clock();
	outfile << euro_spread_price(240) << endl;
	clock_t endTime5 = clock();
	clock_t clockTicksTaken5 = endTime5 - startTime5;
	outfile2 << clockTicksTaken5 / (double)CLOCKS_PER_SEC << endl;

	clock_t startTime6 = clock();
	outfile << euro_spread_price(480) << endl;
	clock_t endTime6 = clock();
	clock_t clockTicksTaken6 = endTime6 - startTime6;
	outfile2 << clockTicksTaken6 / (double)CLOCKS_PER_SEC << endl;
	/*for (int division = 15; division <= 480; division = division*2) {
		cout << euro_spread_price(division) << endl;
		outfile << euro_spread_price(division) << endl;
	}*/
}