#ifndef IC_MYOPIC_H
#define IC_MYOPIC_H
#include "IC_Policy.h"
#include "IC_Random.h"
#include "IC_Support.h"
#include <cmath>
using namespace std;

template <class W>
int Myopic (Policy<W>& P) {
	struct X {	
		static double f (double q, Policy<W> P) {
			vector<double> v = ConditionalExp(q,P);
			return v[0];
		}
		static double goldenSectionSearch(double a, double b, double c, double tau, Policy<W>& P) {
		    double x;
		    double phi = 1.61803398874989; //Golden ratio
			double resphi = 2 - phi; //Golden ratio
		    if (c - b > b - a)
		      x = b + resphi * (c - b);
		    else
		      x = b - resphi * (b - a);
		    if (abs(c - a) < tau * (abs(b) + abs(x))) 
		      return (c + a) / 2; 
		    if (f(x,P) < f(b,P)) {
		      if (c - b > b - a) return goldenSectionSearch(b, x, c, tau, P);
		      else return goldenSectionSearch(a, x, b, tau,P);
		    }
		    else {
		      if (c - b > b - a) return goldenSectionSearch(a, b, x, tau,P);
		      else return goldenSectionSearch(x, b, c, tau,P);
		    }
		}
	};
	while(P.cur<P.T){
		
		P.cur =P.cur+1;
		if (P.cur > P.T - P.L) {
			P.q[P.cur] = 0;
		} 
		else {
			vector<double> v;
			vector<double> w;
			double a = 0;
			double b = 20;
			double c = 40;
			while(true) {
				v = ConditionalExp(b,P);
				w = ConditionalExp(c,P);
				if (v[0] > w[0]) 
					c *= 2;
				else
					break;
			}
			double eps = 0.00001;
			double order = X::goldenSectionSearch(a,b,c,eps,P);
			if (P.is_integer) {
				double x = Unif();
				if (x > order - (int)order) P.q[P.cur] = (int)order;
				else P.q[P.cur] = (int)order + 1;

			} else P.q[P.cur] = order;
		}
		RandomDemand(P);
		P.update();
	}
	return 1;
}


#endif