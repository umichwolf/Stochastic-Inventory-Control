#ifndef IC_DUALBALANCING_H
#define IC_DUALBALANCING_H
#include "IC_Policy.h"
#include "IC_Random.h"
#include "IC_Support.h"
#include <cmath>
using namespace std;

template <class W>
int DualBalancing (Policy<W>& P) {
	while(P.cur<P.T){
		//cout<<"P.cur "<<P.cur<<endl;
		P.cur =P.cur+1;
		if (P.cur > P.T - P.L) {
			P.q[P.cur] = 0;
		} 
		else {
			vector<double> v;
			double min = 0;
			double max = 1;
			while(true) {
				v = ConditionalExp(max,P); // increasing 
				//cout<<v[0]<<' '<<v[1]<<endl;
				if (v[0] - v[1] < 0) 
					max *= 2;
				else
					break;
			}
			//cout<<"updating successfully :"<<P.cur<<endl;
			double eps = 0.0000001;
			double order = (min + max) / 2;
			while (true) {
				v = ConditionalExp(order,P);
				if (abs(v[0] - v[1]) < eps) break;
				else if(v[0] - v[1] > 0) {
					max = order;
					order = (min + order)/2;
				} else {
					min = order;
					order = (max+order)/2;
				}
			}
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