#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>
using namespace std;

double vq1_1 = 0, vq1_2 = 0, pq2_1 = 0, rsq = 0;
double q21_ae, q21_se,q22_se,q22_ae,q31_se,q32_se;
unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);
uniform_real_distribution <double> distribution(0.0, 1.0);
vector <double> european_call; vector <double> dno_call; vector <double> q21_error;
vector <double> q22_error;
double lr_delta_90 = 0, lr_delta_100 = 0, lr_delta_110 = 0;
double lr_gamma_90 = 0, lr_gamma_100 = 0, lr_gamma_110 = 0;
int no_of_simulation = 0;
double *themu, *asian_option;

double max(double a, double b) {
	return (b < a) ? a : b;
}

/*double get_uniform()
{
	return (((double)rand()) / RAND_MAX);
}*/

double get_uniform()
{
	double a = distribution(generator);
	return a;
}

double get_gaussian()
{
	return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}

double get_gaussian_given_u(double u)
{
	return (sqrt(-2.0*log(u))*cos(6.283185307999998*u));
}

double get_anti_gaussian_given_u(double u)
{
	double anti = 1 - u;
	return (sqrt(-2.0*log(anti))*cos(6.283185307999998*anti));
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; };
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double get_mean(vector <double> price, int num) {
	double mean = 0;
	for (int i = 0; i < num; i++) {
		mean = mean + price[i];
	}
	mean = (double)mean / num;
	return mean;
}

double variance(vector <double> price, int num) {
	double mean = 0;
	double variance = 0;
	for (int i = 0; i < num; i++) {
		mean = mean + price[i];
	}
	mean = (double)mean / num;
	for (int i = 0; i < num; i++) {
		variance = variance + pow(price[i] - mean, 2);
	}
	//variance = variance / (num* (num - 1));
	variance = variance / (num);
	return variance;
}

double simq1_2() {
	vector <double> res;
	double v = 0;
	for (int tr = 0; tr < 10000; tr++) {
		/*double x = get_uniform();
		double xt = sqrt(-2.0*log(x)) * cos(6.283185307999998*x);
		while (isfinite(xt) != 1) {
			x = get_uniform(); 
			xt = sqrt(-2.0*log(x)) * cos(6.283185307999998*x);
		}*/
		double temp = 1 / (1 + exp(1 - get_gaussian()));
		res.push_back(temp);
		v += temp;
	}
	
	double vavg = v / 10000; vq1_1 = vavg;
	double var_mc = variance(res,10000);
	return sqrt(var_mc)/sqrt(10000);
}


double simq1_3() {
	vector <double> res1; vector <double> res2; vector <double> w;
	double v1 = 0; double v2 = 0;
	for (int tr = 0; tr < 5000; tr++) {
		double x = get_uniform(); /*double y = 1 - x;
		double xt = sqrt(-2.0*log(x)) * cos(6.283185307999998*x);
		double xt2 = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		while (isfinite(xt) != 1 && isfinite(xt2) != 1) {
			x = get_uniform(); y = 1 - x;
			xt = sqrt(-2.0*log(x)) * cos(6.283185307999998*x);
			xt2 = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		}*/
		double temp = 1.0 / (1 + exp(1 - get_gaussian_given_u(x))); 
		res1.push_back(temp);
		double temp2 = 1.0 / (1 + exp(1 + get_gaussian_given_u(x))); 
		res2.push_back(temp2);
		v1 += temp;
		v2 += temp2;
		w.push_back(0.5*temp + 0.5*temp2);
	}
	double vavg1 = v1 / 5000; double vavg2 = v2 / 5000;
	vq1_2 = 0.5*(vavg1 + vavg2);
	
	double var_mcav = variance(w,w.size());
	return sqrt(var_mcav)/sqrt(5000);
}

double simq2_1() {
	vector <double> res; double allprice = 0;
	int no_of_sim = 10000;
	double ae = 0;
	for (int i = 0; i < no_of_sim; i++) {
		double S = 99, r = 0.03, sigma = 0.6, T = 25, K = 105, B = 90, h = 1.0/252;
		double stay = 1;
		for (int t = 0; t < 24; t++) {
			double R = (r - 0.5*pow(sigma, 2))*h;
			double SD = sigma*sqrt(h);
			S = S* exp(R + SD*get_gaussian());
			if (S <= B) {
				stay = 0;
			}
		}
		//price of each trial
		double one_res = max(0.0, S - K);
		// price of vanilla call
		european_call.push_back(one_res);
		if (stay == 0) {
			one_res = 0;
		}
		// dno call
		res.push_back(one_res);
		dno_call.push_back(one_res);
		ae += (one_res - 4.64765)*(one_res - 4.64765);
		allprice += one_res;
		//actual error of each trial
	}
	pq2_1 = exp(-0.03*25.0/252)*( allprice / (double)no_of_sim);	//avg price
	q21_ae = sqrt(ae / (double) no_of_sim) / sqrt(no_of_sim);
	q21_se = sqrt(variance(res, no_of_sim)/no_of_sim);
	return pq2_1; //price of all trials
}

/*double vanilla_European_call() {
	vector <double> res; double allprice = 0;
	int no_of_sim = 20000;
	for (int i = 0; i < no_of_sim; i++) {
		double S = 99, r = 0.03, sigma = 0.6, T = 25, K = 105, B = 90, h = 1.0 / 252;
		for (int t = 0; t < 24; t++) {
			double R = (r - 0.5*pow(sigma, 2))*h;
			double SD = sigma*sqrt(h);
			S = S* exp(R + SD*get_gaussian());
		}
		//price of each trial
		double one_res = max(0.0, S - K);
		european_call.push_back(one_res);
		allprice += one_res;
	}
	double european_price = exp(-0.03*25.0 / 252)*(allprice / (double)no_of_sim);	//avg price
	return european_price; //price of all trials
}*/

double simq2_2() {
//	Z = X + c(Y - E[Y]), E[Z] = E[X]
	double ymean = pq2_1;
	double allprice = 0;
	int no_of_sim = 10000;
	// make sure pq2_1 is computed before this
	double xmean = get_mean(european_call,no_of_sim);
	double b=0, upper=0, lower=0, ae = 0,testse = 0;
	for (int i = 0; i < no_of_sim; i++) {
		upper += (european_call[i] - xmean)*(dno_call[i] - ymean);
		lower += (european_call[i] - xmean)*(european_call[i] - xmean);
	}
	b = upper / (double) lower;
	double ybar = get_mean(dno_call, no_of_sim);
	vector <double> cvest; double cvres = 0, cvtot = 0;
	for (int i = 0; i < no_of_sim; i++) {
		double temp = dno_call[i] - b*(european_call[i] - xmean);
		allprice += temp;
		cvest.push_back(temp);
		ae += (temp - 4.64765)*(temp - 4.64765);
		//cvres += pow(dno_call[i] - temp, 2);
		//cvtot += pow(dno_call[i] - ybar, 2);
	}
	double est_cv = allprice / (double)no_of_sim;
	//double cvres = 0, cvtot = 0;
	for (int i = 0; i < no_of_sim; i++) {
		//cvres += pow(dno_call[i] - cvest[i], 2);
		//cvtot += pow(dno_call[i] - ybar, 2);
		//cvres += pow( cvest[i]- ybar , 2);
		//cvtot += pow(dno_call[i] - ybar, 2);
		cvres = cvres + pow(european_call[i] - dno_call[i], 2);
		cvtot = cvtot + pow(european_call[i] - xmean, 2);
		//cout << dno_call[i]  << endl; cout << cvest[i] << endl;
		
	}
	q22_se = sqrt(variance(cvest,no_of_sim))/sqrt(no_of_sim);
	q22_ae = sqrt(ae / (double)no_of_sim) / sqrt(no_of_sim);
	//q22_se = sqrt((testse / (double)no_of_sim)) / sqrt(no_of_sim);
	rsq = 1-( cvres / (double)cvtot);
	return est_cv;
}

double vanilla_European_put_delta_pathwise(double K) {
	vector <double> res; double alldelta = 0, lrdelta = 0, lrgammasum =0 ;
	int no_of_sim = 20000;
	for (int i = 0; i < no_of_sim; i++) {
		double S = 100, r = 0.05, sigma = 0.3, T = 0.1, h = T / 30;
		for (int t = 0; t < 30; t++) {
			double R = (r - 0.5*pow(sigma, 2))*h;
			double SD = sigma*sqrt(h);
			S = S* exp(R + SD*get_gaussian());
		}
		//price of each trial
		double payoff = max(0.0, K - S);
		double delta = 0,lrmethod = 0,lrgamma = 0;
		double zeta = (log(S / 100) - (r - 0.5*sigma*sigma)*T) / ( sigma*sqrt(T));
		lrmethod = exp(-r*T)* max(0.0, K - S)*zeta/ (100*sigma*sqrt(T));
		//lrgamma = -exp(-r*T)* max(0.0, K - S)*((zeta*zeta-1)/(100*100* sigma*sqrt(T)*sigma*sqrt(T))- zeta/ (100*100* sigma*sqrt(T)));
		lrgamma = exp(-r*T)* max(0.0, K - S)*((zeta*zeta - zeta*sigma*sqrt(T) -1) / (100 * 100 * sigma* sigma*T));
		if (payoff > 0) {
			delta = -exp(-r*T)*S / 100;	
		}
		alldelta += delta;
		lrdelta += lrmethod; lrgammasum += lrgamma;
	}
	double avg_delta = alldelta / (double)no_of_sim;
	if (K == 90.0) {
		lr_delta_90 = lrdelta / (double)no_of_sim;
		lr_gamma_90 = lrgammasum / (double)no_of_sim;
	}
	else if (K == 100.0) {
		lr_delta_100 = lrdelta / (double)no_of_sim;
		lr_gamma_100 = lrgammasum / (double)no_of_sim;
	}
	else {
		lr_delta_110 = lrdelta / (double)no_of_sim;
		lr_gamma_110 = lrgammasum / (double)no_of_sim;
	}
	//gamma
	/*
	if (K == 90.0) {
		lr_gamma_90 = lrgammasum / (double)no_of_sim;
	}
	else if (K == 100.0) {
		lr_gamma_100 = lrgammasum / (double)no_of_sim;
	}
	else {
		lr_delta_110 = lrgammasum / (double)no_of_sim;
	}*/
	return avg_delta; //price of all trials
}

double option_price_delta_put_black_scholes(const double& S, // spot price
	const double& K, // Strike (exercise) price,
	const double& r,  // interest rate
	const double& sigma,
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double delta = -N(-d1);
	return delta;
}

double option_price_gamma_put_black_scholes(const double& S, // spot price
	const double& K, // Strike (exercise) price,
	const double& r,  // interest rate
	const double& sigma,
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	//double ga = N(d1) / (S*sigma*time_sqrt);
	double ga = exp(-d1*d1/2)/(S*sigma*time_sqrt*sqrt(2* 3.14159265358979323846));
	return ga;
}

//double *find_y(double S0, double K, double r, double sigma, double T, int m)
double *find_y(double S0, double K, double r, double sigma, int m)
{
	double *z = new double[m];
	double *S = new double[m];
	double error = 0.000001;
	//double h = T /(double) m;
	double h = 1.0/252.0;
	double y = 1;
	double sm = S0, indicator = 1;

	for (int i = 0; i < m+1; i++){
		z[i] = 0.0;
	}
	for (int i = 0; i < m+1; i++){
		S[i] = 0.0;
	}
	S[0] = S0;
	while (indicator > error){
		y = y + 0.0000005;
		z[0] = sigma*sqrt(h)*(y + K) / y;
		sm = S0;
		for (int j = 1; j < m+1; j++){
			S[j] = S[j - 1] * exp((r - 0.5*sigma*sigma)*h + sigma*sqrt(h)*z[j - 1]);
			z[j] = z[j - 1] - sigma*sqrt(h)*S[j] / (m*y);
			sm += S[j];
		}
		sm += S[m - 1] * exp((r - 0.5*sigma*sigma)*h + sigma*sqrt(h)*z[m]);
		indicator = max((sm / m - K - y), -(sm / m - K - y));
	}
	return z;
}


double simq3_1() {
	double s0 = 40, r = 0.05, K = 35, sigma = 0.4, h = 1.0 / 252.0;
	int times = 90;
	double asian_call_price = 0;
	vector <double> asians;
	vector <double> path;
	for (int tr = 0; tr < no_of_simulation; tr++) {
		double S = 40; double sum = 40;
		for (int i = 0; i < times; i++) {
			double R = (r - 0.5*pow(sigma, 2))*h;
			double SD = sigma*sqrt(h);
			S = S* exp(R + SD*get_gaussian()); sum += S;
		}
		double asiancall = max(0.0, sum / (double)(times + 1) - K);
		asians.push_back(asiancall);
		asian_call_price += asiancall;
	}
	double asianprice = exp(-r*90.0/252)*asian_call_price / (double) no_of_simulation;
	double se_asian = sqrt(variance(asians, no_of_simulation)/ (double) no_of_simulation);
	q31_se = se_asian;
	return asianprice;
}

double simq3_2(int no_of_sim) {
	double s0 = 40, r = 0.05, K = 35, sigma = 0.4, h = 1.0 / 252.0;
	int times = 90;
	double asian_call_price = 0;
	vector <double> asians;
	
	for (int tr = 0; tr < no_of_sim; tr++) {
		double S = 40; double sum = 40;
		double new_term = 0;
		for (int i = 1; i < times+1; i++) {
			double z =  themu[i]+ get_gaussian();
			new_term = new_term -z * themu[i] + 0.5*pow(themu[i], 2);
			S = S*exp((r - 0.5 * pow(sigma, 2))*h + sigma * sqrt(h)*z);
			sum += S;

		}
		double asiancall = max(0.0, sum / (double)(times + 1) - K)*exp(new_term);
		//	sum_mu*(get_gaussian()+sum_mu)+0.5*sum_mu2);
		asians.push_back(asiancall);
		asian_call_price += asiancall;
	}
	double asianprice = exp(-r*90.0 / 252.0)*asian_call_price / (double) no_of_sim;
	double se_asian = sqrt(variance(asians, no_of_sim) / (double) no_of_sim);
	q32_se = se_asian;
	return asianprice;

}


int main() {
	/*cout << "sdq2   " << simq1_2() << endl;
	cout << "sdq3   "<< simq1_3() << endl;
	cout << "vq2    " << vq1_1 << endl;
	cout << "vq3    " << vq1_2 << endl;
	//cout << "----------------------" << endl;*/
	/*
	cout << "est of Q2-1        "<<  simq2_1() << endl;
	cout <<"SE of Q2-1        "<< q21_se << endl;
	cout << "AE of Q2-1        " << q21_ae << endl;
	cout << "----------------------" << endl;
	cout << "est of Q2-2        " << simq2_2() << endl;
	cout << "SE of Q2-2        " << q22_se << endl;
	cout << "AE of Q2-2        " << q22_ae << endl;
	cout << "R^2        " << rsq << endl;*/
	
	cout << "K90  pathwise delta       " << vanilla_European_put_delta_pathwise(90.0) << endl;
	cout << "K100 pathwise delta       " << vanilla_European_put_delta_pathwise(100.0) << endl;
	cout << "K110 pathwise delta       " << vanilla_European_put_delta_pathwise(110.0) << endl;
	cout << "K90  lr       delta       " << lr_delta_90 << endl;
	cout << "K100 lr       delta       " << lr_delta_100 << endl;
	cout << "K110 lr       delta       " << lr_delta_110 << endl;
	cout << "K90   BS      delta       " << option_price_delta_put_black_scholes(100,90,0.05,0.3,0.1) << endl;
	cout << "K100  BS      delta       " << option_price_delta_put_black_scholes(100, 100, 0.05, 0.3, 0.1) << endl;
	cout << "K110  BS      delta       " << option_price_delta_put_black_scholes(100, 110, 0.05, 0.3, 0.1) << endl;
	cout << endl;
	cout << "K90  lr       gamma       " << lr_gamma_90 << endl;
	cout << "K100 lr       gamma       " << lr_gamma_100 << endl;
	cout << "K110 lr       gamma       " << lr_gamma_110 << endl;
	cout << "K90   BS      gamma       " << option_price_gamma_put_black_scholes(100, 90, 0.05, 0.3, 0.1) << endl;
	cout << "K100  BS      gamma       " << option_price_gamma_put_black_scholes(100, 100, 0.05, 0.3, 0.1) << endl;
	cout << "K110  BS      gamma       " << option_price_gamma_put_black_scholes(100, 110, 0.05, 0.3, 0.1) << endl; 


	/*ofstream output_file("q31");
	for (int i = 100; i <= 10000; i += 50) {
		no_of_simulation = i;
		output_file << simq3_1() << endl;
	}
	cout << q31_se << endl;*/
	
	//cout << mu << endl;
	/*no_of_simulation = 10000;
	cout << simq3_1() << endl;
	cout << q31_se << endl;
	themu = find_y(40.0, 35.0, 0.05, 0.4, 90);*/
	/*cout << simq3_2(10000) << endl;
	cout << q32_se << endl;*/
	/*
	for (int i = 0; i <= 100; i++) {
		cout << simq3_2(10000) << endl;
		cout << q32_se << endl;
		cout << "...." << endl;
	}*/
	/*
	ofstream output_file2("q32");
	for (int i = 100; i <= 10000; i += 50) {
		output_file2 << simq3_2(i) << endl;
	}*/
	/*no_of_simulation = 10000;
	cout << simq3_1() << endl;
	cout << q31_se << endl;
	themu = find_y(40.0, 35.0, 0.05, 0.4, 90); 
	cout << simq3_2(10000) << endl;
	cout << q32_se << endl;*/
}