#ifndef IC_PARAMETERDATA_H
#define IC_PARAMETERDATA_H
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

//ParameterData is a class of given parameter for our inventory control problem
class ParameterData{
public:
	//The following are the parameters that we are going to specify before doing the optimization problem
	//Planning Horizon
	int T;
	//Leading Time
	int L;
	//Ordering Cost
	vector<double> c;
	//Holding Cost
	vector<double> h;
	//Backlogging Cost
	vector<double> p;
	//Discount rate
	double alpha;


	//Constructors
	//Default Constructor
	ParameterData();
	//Constructor assigning constant cost coefficient
	ParameterData(int vT,int vL,double vc,double vh,double vp,double valpha);
	//general constructor
	ParameterData(int vT,int vL,vector<double> vc,vector<double> vh,vector<double> vp,double valpha);

	//Ouput the parameters
	int output();
	//Output the parameters at time t
	int output_t(int t);
};

ParameterData::ParameterData(){
	T=0;
	L=0;
	alpha=1;
	
	c.push_back(0);
	h.push_back(0);
	p.push_back(0);
}

ParameterData::ParameterData(int vT,int vL,double vc,double vh,double vp,double valpha){
	int i;
	T=vT;
	L=vL;
	alpha=valpha;

	c.push_back(0);
	h.push_back(0);
	p.push_back(0);

	for(i=1;i<=T;i++){
		c.push_back(vc);
		h.push_back(vh);
		p.push_back(vp);
	}
}

ParameterData::ParameterData(int vT,int vL,vector<double> vc,vector<double> vh,vector<double> vp,double valpha){
	int i;
	T=vT;
	L=vL;
	alpha=valpha;
	
	c.push_back(0);
	h.push_back(0);
	p.push_back(0);

	for(i=1;i<=T;i++){
		c.push_back(vc[i-1]);
		h.push_back(vh[i-1]);
		p.push_back(vp[i-1]);
	}
}

int ParameterData::output(){
	int i;
	cout.precision(7);
	cout<<"Planning Horizon is: "<<T<<endl;
	cout<<"Leading Time is: "<<L<<endl;
	cout<<"Discount Rate is: "<<alpha<<endl;

	cout<<setw(6)<<"Period"<<setw(15)<<"Ordering Cost"<<setw(14)<<"Holding Cost"<<setw(18)<<"Backlogging Cost"<<endl;
	for(i=0;i<=T;i++)
		cout<<setw(6)<<i<<setw(15)<<c[i]<<setw(14)<<h[i]<<setw(18)<<p[i]<<endl;
	cout<<endl;
	return 1;
}

#endif