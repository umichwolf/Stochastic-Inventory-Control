#include "IC_Policy.h"
#include "IC_Random.h"
#include <cmath>     
#include <map>   // std::abs
using namespace std;

int T = 2;
const int n = 3;
double qs[n] = {0,1,2};
double ds[n] = {0,1,2};
double p[n] = {0.1, 0.7, 0.2};

struct State {
	int t,x;
 	bool operator<( const State& other) const {
          if ( t == other.t ) {
          	return x < other.x;
          } else return t < other.t;
      }
};

map <State, double> bag;

double g(int x, int q, int d) { // local cost function
	return q + (x+q-d)*(x+q-d);
}

double dot(double* v, double* w) { // helper function, dot product
	double sum = 0;
	for (int i = 0; i<n; i++) {
		sum += v[i] * w[i];
	}
	return sum;
}

int next_x(int x, int q, int d) { // transition function of x
	return max(0,x+q-d);
}

double J(int t, int x) {
	if (t==T) {
		double q0 = -1;
		double cost0 = 10000000;
		for (int i = 0; i<n; i++) {
			double v[n];
			for (int j = 0; j<n; j++) {
				v[j] = g(x,qs[i],ds[j]);
			}
			double cost = dot(p, v);
			if(cost < cost0) {
				q0 = qs[i];
				cost0 = cost;
			}
		}
		State s = {t,x};
		bag[s] = cost0;
		return cost0;
	} else {
		double q0 = -1;
		double cost0 = 10000000;
		for (int i = 0; i<n; i++) {
			double v[n];
			for (int j = 0; j<n; j++) {
				State s = {t+1,next_x(x,qs[i],ds[j])};
				if (!bag[s]) bag[s] = J(t+1, next_x(x,qs[i],ds[j]));
				v[j] = g(x,qs[i],ds[j]) + bag[s];
				}
			double cost = dot(p,v);
			if(cost < cost0) {
				q0 = qs[i];
				cost0 = cost;
			}
		}
		return cost0;
	}
}
double mu(int x, int t) { // DP policy
	// if (t == T) return
	return 0;
}

int main () {
	cout << J(0,0) << endl;
	return 0;
}
