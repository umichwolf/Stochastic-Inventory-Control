#ifndef IC_POLICY_H
#define IC_POLICY_H
#include"IC_ParameterData.h"
#include"IC_Support.h"
#include<fstream>
#include<algorithm>
using namespace std;

const int ALG=3;
const double ALPHA=1;
const int X0=0;
const int IS_INTEGER=0;
const int ITEMNUM=30;
const int DIST=0;
const vector<double> RD_P(1,0);

//Policy is a class saving result and sample path for different algorithm approaches
template<class W> class Policy;
//The three are the core algorithm to compute the inventory control problem
template<class W> int DualBalancing (Policy<W>& P);
template<class W> int Myopic (Policy<W>& P);
template<class W> int DP (Policy<W>& P);

//The following are the friend function in IC_Random.h
//This function return all the different conditional expectations based on different algorithms and distributions
template<class W> vector<double> ConditionalExp(double q,Policy<W> &policy);
//The function generate a sample of demand based on the distributions
template<class W> int RandomDemand(Policy<W> &policy);
//These returns conditional expectation depends on different distribution
template<class W> vector<double> BinomialCE(double q,Policy<W> &policy);
template<class W> vector<double> IncreBinomialCE(double q,Policy<W> &policy);
template<class W> vector<double> NormalCE(double q,Policy<W> &policy);
template<class W> vector<double> RA1NormalCE(double q, Policy<W> &P);
template<class W> vector<double> RAnNormalCE(double q, Policy<W> &P);

//This class collect datas
template<class W> class Samples;

//Policy is a class saving result and sample path for different algorithm approaches
template <class W>
class Policy:public ParameterData{
public:
	//The following are the datas we are going to generate while doing the optimization problem

	//alg include what algorithm we used to compute the policy
	//1 for DP, 2 for myopic, 3 for dual balancing
	int alg;
	//dist is what kind of random process distribution we used to generat the demand
	int dist;
	//rd_p include all the parameter for the distribution of the random process 
	vector<double> rd_p;
	//whether W is integer type
	int is_integer;

	//Current time period that we already know the inventory position, demand and ordering.
	int cur;
	//Cost Vector
	vector<double> C;
	//Holding Cost Vector
	vector<double> Ch;
	//Backlogging Cost Vector
	vector<double> Cp;
	//Inventory Position
	vector<W> x;
	//Net Inventory
	vector<W> NI;
	//Demand
	vector<W> d;
	//Order
	vector<W> q;
	//Accumulative Cost Vector
	vector<double> C_accum;
	//Accumulative Holding Cost Vector
	vector<double> Ch_accum;
	//Accumulative Backlogging Cost Vector
	vector<double> Cp_accum;
	//Accumulative Demand
	vector<W> d_accum;
	//Accumulative Order
	vector<W> q_accum;

	//Total Cost,Demand and Order
	double total_C;
	double total_Ch;
	double total_Cp;
	W total_d;
	W total_q;
	
	//Constructors
	//Default Constructor
	Policy();
	//Constructor assigning constant cost coefficient
	Policy(int vT,int vL,double vc,double vh,double vp,int valg=ALG,double valpha=ALPHA,W vx0=(W)(X0),int vis_integer=IS_INTEGER,int vdist=DIST,vector<double> vrd_p=RD_P);
	//general constructor
	Policy(int vT,int vL,vector<double> vc,vector<double> vh,vector<double> vp,int valg=ALG,double valpha=ALPHA,W vx0=(W)(X0),int vis_integer=IS_INTEGER,int vdist=DIST,vector<double> vrd_p=RD_P);


	//Out put the sample path or the state set
	int output();
	//Output on screen
	int soutput();
	//Output to file
	int foutput();

	//This function update the sample path or the state after we get every demand and ordering in every 
	int update();
	
	//The following are friend functions or classes so we will be able to access the core data in the Policy class
	//ConditionalExp is a friend function of class Policy
	friend vector<double> ConditionalExp<>(double q,Policy<W> &policy);
	//RandomDemand is a friend function of class Policy
	friend int RandomDemand<>(Policy<W> &policy);
	//DualBalancing/DP.Myopic are friend functions of class Policy
	friend int DualBalancing<>(Policy<W> &policy);
	friend int DP<>(Policy<W> &policy);
	friend int Myopic<>(Policy<W> &policy);

	//Different conditional expectatoins. 
	friend vector<double> BinomialCE<>(double q,Policy<W> &policy);
	friend vector<double> IncreBinomialCE<>(double q,Policy<W> &policy);
	friend vector<double> NormalCE<>(double q,Policy<W> &policy);
	friend vector<double> RA1NormalCE<>(double q, Policy<W> &P);
	friend vector<double> RAnNormalCE<>(double q, Policy<W> &P);

	
	//Samples class is a friend class of Policy
	friend class Samples<W>;
};

template <class W>
Policy<W>::Policy():ParameterData(){
	alg=0;
	dist=0;
	rd_p.push_back(0);
	is_integer=0;

	cur=0;
	total_C=0;
	total_Ch=0;
	total_Cp=0;
	total_d=0;
	total_q=0;

	C.push_back(0);
	Ch.push_back(0);
	Cp.push_back(0);
	x.push_back(0);
	NI.push_back(0);
	d.push_back(0);
	q.push_back(0);

	C_accum.push_back(0);
	Ch_accum.push_back(0);
	Cp_accum.push_back(0);
	d_accum.push_back(0);
	q_accum.push_back(0);
}

template <class W>
Policy<W>::Policy(int vT,int vL,double vc,double vh,double vp,int valg,double valpha,W vx0,int vis_integer,int vdist,vector<double> vrd_p):ParameterData(vT,vL,vc,vh,vp,valpha){
	alg=valg;
	dist=vdist;
	rd_p=vrd_p;
	is_integer=vis_integer;

	cur=0;
	total_C=0;
	total_Ch=0;
	total_Cp=0;
	total_d=0;
	total_q=0;

	C.assign(vT+1,0);
	Ch.assign(vT+1,0);
	Cp.assign(vT+1,0);
	C_accum.assign(vT+1,0);
	Ch_accum.assign(vT+1,0);
	Cp_accum.assign(vT+1,0);

	x.assign(vT+1,0);
	d.assign(vT+1,0);
	q.assign(vT+1,0);
	NI.assign(vT+1,0);
	d_accum.assign(vT+1,0);
	q_accum.assign(vT+1,0);

	x[0]=vx0;
	NI[0]=vx0;
}

template <class W>
Policy<W>::Policy(int vT,int vL,vector<double> vc,vector<double> vh,vector<double> vp,int valg,double valpha,W vx0,int vis_integer,int vdist,vector<double> vrd_p):ParameterData(vT,vL,vc,vh,vp,valpha){
	alg=valg;
	dist=vdist;
	rd_p=vrd_p;
	is_integer=vis_integer;

	cur=0;
	total_C=0;
	total_Ch=0;
	total_Cp=0;
	total_d=0;
	total_q=0;

	C.assign(vT+1,0);
	Ch.assign(vT+1,0);
	Cp.assign(vT+1,0);
	C_accum.assign(vT+1,0);
	Ch_accum.assign(vT+1,0);
	Cp_accum.assign(vT+1,0);

	x.assign(vT+1,0);
	d.assign(vT+1,0);
	q.assign(vT+1,0);
	NI.assign(vT+1,0);
	d_accum.assign(vT+1,0);
	q_accum.assign(vT+1,0);

	x[0]=vx0;
	NI[0]=vx0;
}


template <class W>
int Policy<W>::output(){
	int temp;
	cout<<"Type 1 if you want output to screen, 0 otherwise."<<endl;
	cin>>temp;
	if(temp==1)
		soutput();
	temp=0;
	cout<<"Type 1 if you want output to a file, 0 otherwise."<<endl;
	cin>>temp;
	if(temp==1)
		foutput();
	return 1;
}


template <class W>
int Policy<W>::soutput(){
	int i;
	int rd_p_size=rd_p.size();
	cout.precision(4);
	cout<<"Planning Horizon is: "<<T<<endl;
	cout<<"The Current Step is at time: "<<cur<<endl;
	cout<<"Leading Time is: "<<L<<endl;
	cout<<"Discount Rate is: "<<alpha<<endl;
	cout<<"Initial Inventory is: "<<x[0]<<endl;
	cout<<"This is a integer typed problem or not: "<<is_integer<<endl;
	cout<<"The Algorithm is: "<<alg<<endl;
	cout<<"The Distribution Type is: "<<dist<<endl;
	cout<<"The Distribution Parameters are: ";
	for(i=0;i<rd_p_size;i++)
		cout<<setw(10)<<rd_p[i];
	cout<<endl;
	cout<<"Total Cost is: "<<total_C<<endl;
	cout<<"Total Demand is: "<<total_d<<endl;
	cout<<"Total Order is: "<<total_q<<endl;
	cout<<"Total Holding Cost is: "<<total_Ch<<endl;
	cout<<"Total Backlogging Cost is: "<<total_Cp<<endl;
	cout<<endl;

	i=0;
	getchar();
	while(1){
		cout<<setw(6)<<"Period"<<setw(11)<<"Inv Pos"<<setw(11)<<"Net Inv"<<setw(11)<<"Cost"<<setw(11)<<"Demand"<<setw(11)<<"Ordering"<<setw(11)<<"H Cost"<<setw(11)<<"B Cost";
		cout<<setw(11)<<"Cost Ac"<<setw(11)<<"Demand Ac"<<setw(11)<<"Order Ac"<<setw(11)<<"H Cost Ac"<<setw(11)<<"B Cost Ac"<<endl;
		for(;i<=T;i++){
			cout<<setw(6)<<i<<setw(11)<<x[i]<<setw(11)<<NI[i]<<setw(11)<<C[i]<<setw(11)<<d[i]<<setw(11)<<q[i]<<setw(11)<<Ch[i]<<setw(11)<<Cp[i];
			cout<<setw(11)<<C_accum[i]<<setw(11)<<d_accum[i]<<setw(11)<<q_accum[i]<<setw(11)<<Ch_accum[i]<<setw(11)<<Cp_accum[i]<<endl;
			if(i%ITEMNUM==0&&i!=0)
				break;
		}
		cout<<endl;
		if(i>=T)
			break;
		else{
			cout<<"Press Enter to continue."<<endl;
			getchar();
			i=i+1;
		}
	}
	return 1;
}

template <class W>
int Policy<W>::foutput(){
	char filename[200];
	cout<<"Please type a file name for output. (.txt)"<<endl;
	while(1){
		cin>>filename;
		ofstream fout;
		fout.open(filename,ios::out);
		if(fout){
			fout.seekp(0);
			cout<<"File is openned successfully, outputting data..."<<endl;
					
			int i;
			int rd_p_size=rd_p.size();
			fout.precision(4);
			fout<<"Planning Horizon is: "<<T<<endl;
			fout<<"The Current Step is at time: "<<cur<<endl;
			fout<<"Leading Time is: "<<L<<endl;
			fout<<"Discount Rate is: "<<alpha<<endl;
			fout<<"Initial Inventory is: "<<x[0]<<endl;
			fout<<"This is a integer typed problem or not: "<<is_integer<<endl;
			fout<<"The Algorithm is: "<<alg<<endl;
			fout<<"The Distribution Type is: "<<dist<<endl;
			fout<<"The Distribution Parameters are: ";
			for(i=0;i<rd_p_size;i++)
				fout<<setw(10)<<rd_p[i];
			fout<<endl;
			fout<<"Total Cost is: "<<total_C<<endl;
			fout<<"Total Demand is: "<<total_d<<endl;
			fout<<"Total Order is: "<<total_q<<endl;
			fout<<"Total Holding Cost is: "<<total_Ch<<endl;
			fout<<"Total Backlogging Cost is: "<<total_Cp<<endl;
			fout<<endl;

			
			fout<<setw(6)<<"Period"<<setw(11)<<"Inv Pos"<<setw(11)<<"Net Inv"<<setw(11)<<"Cost"<<setw(11)<<"Demand"<<setw(11)<<"Ordering"<<setw(11)<<"H Cost"<<setw(11)<<"B Cost";
			fout<<setw(11)<<"Cost Ac"<<setw(11)<<"Demand Ac"<<setw(11)<<"Order Ac"<<setw(11)<<"H Cost Ac"<<setw(11)<<"B Cost Ac"<<endl;
			for(i=0;i<=T;i++){
				fout<<setw(6)<<i<<setw(11)<<x[i]<<setw(11)<<NI[i]<<setw(11)<<C[i]<<setw(11)<<d[i]<<setw(11)<<q[i]<<setw(11)<<Ch[i]<<setw(11)<<Cp[i];
				fout<<setw(11)<<C_accum[i]<<setw(11)<<d_accum[i]<<setw(11)<<q_accum[i]<<setw(11)<<Ch_accum[i]<<setw(11)<<Cp_accum[i]<<endl;
			}
			fout<<endl;

			fout.close();
			cout<<"Data output finished."<<endl<<endl;
			break;
		}
	}
	return 1;
}

template<class W>
int Policy<W>::update(){
	if(cur<=T){
		if(cur<L+1)
			NI[cur]=NI[cur-1]-d[cur];
		if(cur>L)
			NI[cur]=NI[cur-1]-d[cur]+q[cur-L];
		

		x[cur]=x[cur-1]+q[cur]-d[cur];

		Ch[cur]=h[cur]*max((double)NI[cur],0.0);
		Cp[cur]=p[cur]*max(-(double)NI[cur],0.0);
		C[cur]=Ch[cur]+Cp[cur]+c[cur]*q[cur];
		
		C_accum[cur]=C_accum[cur-1]+C[cur];
		Ch_accum[cur]=Ch_accum[cur-1]+Ch[cur];
		Cp_accum[cur]=Cp_accum[cur-1]+Cp[cur];
		d_accum[cur]=d_accum[cur-1]+d[cur];
		q_accum[cur]=q_accum[cur-1]+q[cur];
		
		total_C=C_accum[cur];
		total_Ch=Ch_accum[cur];
		total_Cp=Cp_accum[cur];
		total_d=d_accum[cur];
		total_q=q_accum[cur];
	}	
	return 1;
}



#endif