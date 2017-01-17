#ifndef IC_STATISTICS_H
#define IC_STATISTICS_H
#include"IC_Policy.h"
#include"IC_DP.h"
#include"IC_DualBalancing.h"
#include"IC_Myopic.h"
#include"IC_Support.h"
using namespace std;

double duration=0;

//Data Collection Function
int DataCollector();

template<class W>
class Samples{
public:
	//This is the policy class that on which we are going to run the algorithm many times to collect datas
	Policy<W> policy;
	//N is how many times we wrong the algorithms
	int N;

	//The following are parameters determined at very beginning
	//Planning Horizon
	int T;
	//Leading Time
	int L;
	//Discount rate
	double alpha;

	//The following are the datas we are going to generate while doing the optimization problem MANY TIMES

	//alg include what algorithm we used to compute the policy
	//1 for DP, 2 for myopic, 3 for dual balancing
	int alg;
	//dist is what kind of random process distribution we used to generat the demand
	int dist;
	//rd_p include all the parameter for the distribution of the random process 
	vector<double> rd_p;
	//whether W is integer type
	int is_integer;

	//Cost Vector
	vector<double> C;
	//Holding Cost Vector
	vector<double> Ch;
	//Backlogging Cost Vector
	vector<double> Cp;
	//Inventory Position
	vector<double> x;
	//Net Inventory
	vector<double> NI;
	//Demand
	vector<double> d;
	//Order
	vector<double> q;
	//Accumulative Cost Vector
	vector<double> C_accum;
	//Accumulative Holding Cost Vector
	vector<double> Ch_accum;
	//Accumulative Backlogging Cost Vector
	vector<double> Cp_accum;
	//Accumulative Demand
	vector<double> d_accum;
	//Accumulative Order
	vector<double> q_accum;

	//Total Cost,Demand and Order
	double total_C;
	double total_Ch;
	double total_Cp;
	double total_d;
	double total_q;

	//Compute the policy sample once and collect one more data
	int Update();

	//Generate the average of the data given N and policy
	int AverageGenerator();

	//Synchronize the data of the policy and Sample. 
	int Synchronize();


	//Constructors
	//The default constructor
	Samples();
	//Constructor assigning constant cost coefficient
	Samples(int vN,int vT,int vL,double vc,double vh,double vp,int valg=ALG,double valpha=ALPHA,W vx0=(W)(X0),int vis_integer=IS_INTEGER,int vdist=DIST,vector<double> vrd_p=RD_P);
	//general constructor
	Samples(int vN,int vT,int vL,vector<double> vc,vector<double> vh,vector<double> vp,int valg=ALG,double valpha=ALPHA,W vx0=(W)(X0),int vis_integer=IS_INTEGER,int vdist=DIST,vector<double> vrd_p=RD_P);
	//constructor using a Policy class
	Samples(int vN,Policy<W> vpolicy);
	

	//Out put the samples data
	int output();
	//Output on screen
	int soutput();
	//Output to file
	int foutput();
	//Automatic output to file
	int autofoutput(int j);
	//Output a vector of main results
	vector<double> outputmain();

	//DataCollector is a friend function
	friend int DataCollector();

};

template<class W>
int Samples<W>::Update(){
	int i;
	if(alg==1)
		DP(policy);
	if(alg==2)
		Myopic(policy);
	if(alg==3)
		DualBalancing(policy);
	
	for(i=1;i<=T;i++){
		C[i]=C[i]+policy.C[i]/N;
		Cp[i]=Cp[i]+policy.Cp[i]/N;
		Ch[i]=Ch[i]+policy.Ch[i]/N;
		x[i]=x[i]+((double)policy.x[i])/N;
		d[i]=d[i]+((double)policy.d[i])/N;
		q[i]=q[i]+((double)policy.q[i])/N;
		NI[i]=NI[i]+abs(((double)policy.NI[i])/N);
		C_accum[i]=C_accum[i]+policy.C_accum[i]/N;
		Ch_accum[i]=Ch_accum[i]+policy.Ch_accum[i]/N;
		Cp_accum[i]=Cp_accum[i]+policy.Cp_accum[i]/N;
		d_accum[i]=d_accum[i]+((double)policy.d_accum[i])/N;
		q_accum[i]=q_accum[i]+((double)policy.q_accum[i])/N;
	}
	return 1;
}

template<class W>
int Samples<W>::AverageGenerator(){
	int i;
	double start,finish;
	start=clock();
	for(i=1;i<=N;i++){
		Update();
		policy.cur=0;
		cout<<'\r'<<"Computing the "<<setw(5)<<i<<"th simulation for T="<<setw(5)<<T;
	}
	finish=clock();
	duration=(double)(finish-start)/CLOCKS_PER_SEC;
	total_C=C_accum[T];
	total_Ch=Ch_accum[T];
	total_Cp=Cp_accum[T];
	total_d=d_accum[T];
	total_q=q_accum[T];
	return 1;
}


template<class W>
int Samples<W>::Synchronize(){
	T=policy.T;
	L=policy.L;
	alpha=policy.alpha;

	alg=policy.alg;
	dist=policy.dist;
	rd_p=policy.rd_p;
	is_integer=policy.is_integer;

	total_C=policy.total_C;
	total_Ch=policy.total_Ch;
	total_Cp=policy.total_Cp;
	total_d=policy.total_d;
	total_q=policy.total_q;

	C.assign(T+1,0);
	Ch.assign(T+1,0);
	Cp.assign(T+1,0);
	C_accum.assign(T+1,0);
	Ch_accum.assign(T+1,0);
	Cp_accum.assign(T+1,0);

	x.assign(T+1,0);
	d.assign(T+1,0);
	q.assign(T+1,0);
	NI.assign(T+1,0);
	d_accum.assign(T+1,0);
	q_accum.assign(T+1,0);
	x[0]=policy.x[0];
	NI[0]=policy.NI[0];
	return 1;
}


template<class W>
Samples<W>::Samples(){
	Policy<W> vpolicy;
	policy=vpolicy;
	N=0;
	Synchronize();
}

template<class W>
Samples<W>::Samples(int vN,int vT,int vL,double vc,double vh,double vp,int valg,double valpha,W vx0,int vis_integer,int vdist,vector<double> vrd_p){
	Policy<W> vpolicy(vT,vL,vc,vh,vp,valg,valpha,vx0,vis_integer,vdist,vrd_p);
	policy=vpolicy;
	N=vN;
	Synchronize();
}

template<class W>
Samples<W>::Samples(int vN,int vT,int vL,vector<double> vc,vector<double> vh,vector<double> vp,int valg,double valpha,W vx0,int vis_integer,int vdist,vector<double> vrd_p){
	Policy<W> vpolicy(vT,vL,vc,vh,vp,valg,valpha,vx0,vis_integer,vdist,vrd_p);
	policy=vpolicy;
	N=vN;
	Synchronize();
}


template<class W>
Samples<W>::Samples(int vN, Policy<W> vpolicy){
	policy=vpolicy;
	N=vN;
	Synchronize();
}



template <class W>
int Samples<W>::output(){
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
int Samples<W>::soutput(){
	int i;
	int rd_p_size=rd_p.size();
	cout.precision(4);
	cout<<"The data is average of running the code "<<N<<" many times."<<endl;
	cout<<"Planning Horizon is: "<<T<<endl;
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
int Samples<W>::foutput(){
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
			fout<<"The data is average of running the code "<<N<<" many times."<<endl;
			fout<<"Planning Horizon is: "<<T<<endl;
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

template <class W>
int Samples<W>::autofoutput(int j){
	char filename[20]={'o','u','t','p','u','t'};
	char c=j/10+48;
	filename[6]=c;
	c=j%10+48;
	filename[7]=c;
	filename[8]='.';
	filename[9]='t';
	filename[10]='x';
	filename[11]='t';
	filename[12]='\0';

	while(1){
		ofstream fout;
		fout.open(filename,ios::out);
		if(fout){
			fout.seekp(0);
			cout<<'\r'<<filename<<" is openned successfully, outputting data......"<<endl;
			
			int i;
			int rd_p_size=rd_p.size();
			fout.precision(4);
			cout<<"The time it takes to compute per each caseis "<<duration/N<<" secs"<<endl;
			fout<<"The time it takes to compute per each caseis "<<duration/N<<" secs"<<endl;
			fout<<"The data is average of running the code "<<N<<" many times."<<endl;
			fout<<"Planning Horizon is: "<<T<<endl;
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
/*
			fout<<"Period"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<i;
			fout<<endl<<endl;
			fout<<"Inventory Position"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<x[i];
			fout<<endl<<endl;
			fout<<"Net Inventory"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<NI[i];
			fout<<endl<<endl;
			fout<<"Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<C[i];
			fout<<endl<<endl;
			fout<<"Demand"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<d[i];
			fout<<endl<<endl;
			fout<<"Ordering"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<q[i];
			fout<<endl<<endl;
			fout<<"Holding Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<Ch[i];
			fout<<endl<<endl;
			fout<<"Backlogging Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<Cp[i];
			fout<<endl<<endl;
			fout<<"Accumulative Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<C_accum[i];
			fout<<endl<<endl;
			fout<<"Accumulative Demand"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<d_accum[i];
			fout<<endl<<endl;
			fout<<"Accumulative Ordering"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<q_accum[i];
			fout<<endl<<endl;
			fout<<"Accumulative Holding Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<Ch_accum[i];
			fout<<endl<<endl;
			fout<<"Accumulative Backlogging Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<setw(11)<<Cp_accum[i];
			fout<<endl<<endl;
*/
			fout<<"Period"<<endl;
			for(i=0;i<=T;i++)
				fout<<i<<endl;
			fout<<endl<<endl;
			fout<<"InventoryPosition"<<endl;
			for(i=0;i<=T;i++)
				fout<<x[i]<<endl;
			fout<<endl<<endl;
			fout<<"NetInventory"<<endl;
			for(i=0;i<=T;i++)
				fout<<NI[i]<<endl;
			fout<<endl<<endl;
			fout<<"Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<C[i]<<endl;
			fout<<endl<<endl;
			fout<<"Demand"<<endl;
			for(i=0;i<=T;i++)
				fout<<d[i]<<endl;
			fout<<endl<<endl;
			fout<<"Ordering"<<endl;
			for(i=0;i<=T;i++)
				fout<<q[i]<<endl;
			fout<<endl<<endl;
			fout<<"HoldingCost"<<endl;
			for(i=0;i<=T;i++)
				fout<<Ch[i]<<endl;
			fout<<endl<<endl;
			fout<<"BackloggingCost"<<endl;
			for(i=0;i<=T;i++)
				fout<<Cp[i]<<endl;
			fout<<endl<<endl;
			fout<<"AccumulativeCost"<<endl;
			for(i=0;i<=T;i++)
				fout<<C_accum[i]<<endl;
			fout<<endl<<endl;
			fout<<"AccumulativeDemand"<<endl;
			for(i=0;i<=T;i++)
				fout<<d_accum[i]<<endl;
			fout<<endl<<endl;
			fout<<"AccumulativeOrdering"<<endl;
			for(i=0;i<=T;i++)
				fout<<q_accum[i]<<endl;
			fout<<endl<<endl;
			fout<<"AccumulativeHolding Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<Ch_accum[i]<<endl;
			fout<<endl<<endl;
			fout<<"AccumulativeBacklogging Cost"<<endl;
			for(i=0;i<=T;i++)
				fout<<Cp_accum[i]<<endl;
			fout<<endl<<endl;

			fout.close();
			cout<<filename<<" data output finished."<<endl<<endl;
			break;
		}
	}
	return 1;
}

template<class W>
vector<double> Samples<W>::outputmain(){
	vector<double> temp;
	temp.push_back(total_C);
	temp.push_back(total_Ch);
	temp.push_back(total_Cp);
	temp.push_back((double)(total_d));
	temp.push_back((double)(total_q));
	return temp;
}



int DataCollector(){
	char inputfile[200];
	DC_Exp();

	cout<<"Now, please type your input file name. (.txt)"<<endl;
	cin>>inputfile;

	ifstream fin(inputfile,ios::ate);
	if(fin){
		cout<<"The inputfile successfully loaded."<<endl;
		fin.seekg(0);

		int datanum;

		int datatype;
		int vN,vT,vL,vis_integer,valg;
		vector<int> vTv;
		int vTv_l;
		double valpha,vx0;
		int vdist;
		vector<double> vrd_p;
		int vrd_p_l;
		int constcoef;
		double vc,vh,vp;
		vector<double> vcv,vhv,vpv;
		
		fin>>datanum;
		cout<<"Data parameter loading started."<<endl<<endl;

		int i,j,temp;
		double tempd;
		for(j=1;j<=datanum;j++){
			//cout<<"the datanum is "<<j<<endl;
			cout<<"Now, loading the "<<j<<"th data parameter sets."<<endl;
			fin>>datatype;
			fin>>vN;
			if(datatype==1)
				fin>>vT;
			if(datatype==2){
				vTv.clear();
				fin>>vTv_l;
				for(i=1;i<=vTv_l;i++){
					fin>>temp;
					vTv.push_back(temp);
				}
			}
			
			
			fin>>vL>>valpha>>vx0>>vis_integer>>valg>>vdist>>vrd_p_l;
			//cout<<vx0<<endl;
			vrd_p.clear();
			for(i=1;i<=vrd_p_l;i++){
				fin>>tempd;
				vrd_p.push_back(tempd);
			}
			fin>>constcoef;
			if(constcoef==1){
				fin>>vc>>vh>>vp;
				//cout<<"vc,vh,vp"<<vc<<' '<<vh<<' '<<vp<<endl;
			}
			else if(constcoef==0){
				vcv.clear();
				vhv.clear();
				vpv.clear();
				for(i=1;i<=vT;i++){
					fin>>tempd;
					vcv.push_back(tempd);
				}
				for(i=1;i<=vT;i++){
					fin>>tempd;
					vhv.push_back(tempd);
				}
				for(i=1;i<=vT;i++){
					fin>>tempd;
					vpv.push_back(tempd);
				}
			}
			cout<<j<<"th data parameter sets loading finished, start computation."<<endl;

			if(datatype==1){
				if(vis_integer==0){
					if(constcoef==1){
						Samples<double> sample(vN,vT,vL,vc,vh,vp,valg,valpha,vx0,vis_integer,vdist,vrd_p);
						sample.AverageGenerator();
						sample.autofoutput(j);
					}
					else if(constcoef==0){
						Samples<double> sample(vN,vT,vL,vcv,vhv,vpv,valg,valpha,vx0,vis_integer,vdist,vrd_p);
						sample.AverageGenerator();
						sample.autofoutput(j);
					}
				}
				if(vis_integer==1){
					if(constcoef==1){
						Samples<int> sample(vN,vT,vL,vc,vh,vp,valg,valpha,(int)(vx0),vis_integer,vdist,vrd_p);
						sample.AverageGenerator();
						sample.autofoutput(j);
					}
					else if(constcoef==0){
						Samples<int> sample(vN,vT,vL,vcv,vhv,vpv,valg,valpha,(int)(vx0),vis_integer,vdist,vrd_p);
						sample.AverageGenerator();
						sample.autofoutput(j);
					}
				}
			}//end of data type 1

			else if(datatype==2){
				char filename[20]={'o','u','t','p','u','t'};
				char c=j/10+48;
				filename[6]=c;
				c=j%10+48;
				filename[7]=c;
				filename[8]='.';
				filename[9]='t';
				filename[10]='x';
				filename[11]='t';
				filename[12]='\0';
				
				while(1){
					ofstream fout;
					fout.open(filename,ios::out);
					if(fout){
						fout.seekp(0);
						cout<<filename<<" is openned successfully, computing......"<<endl;

						int vrd_p_size=vrd_p.size();
						fout.precision(4);
						fout<<"The data is average of running the code "<<vN<<" many times."<<endl;
						fout<<"Leading Time is: "<<vL<<endl;
						fout<<"Discount Rate is: "<<valpha<<endl;
						fout<<"Initial Inventory is: "<<vx0<<endl;
						fout<<"This is a integer typed problem or not: "<<vis_integer<<endl;
						fout<<"The Algorithm is: "<<valg<<endl;
						fout<<"The Distribution Type is: "<<vdist<<endl;
						fout<<"The Distribution Parameters are: ";
						for(i=0;i<vrd_p_size;i++)
							fout<<setw(10)<<vrd_p[i];
						fout<<endl;

						int ii;
						vector< vector<double> > results;
						for(ii=0;ii<vTv_l;ii++){
							if(vis_integer==0){
								if(constcoef==1){
									Samples<double> sample(vN,vTv[ii],vL,vc,vh,vp,valg,valpha,vx0,vis_integer,vdist,vrd_p);
									sample.AverageGenerator();
									results.push_back(sample.outputmain());
								}
								else if(constcoef==0){
									Samples<double> sample(vN,vTv[ii],vL,vcv,vhv,vpv,valg,valpha,vx0,vis_integer,vdist,vrd_p);
									sample.AverageGenerator();
									results.push_back(sample.outputmain());
								}
							}
							if(vis_integer==1){
								if(constcoef==1){
									Samples<int> sample(vN,vTv[ii],vL,vc,vh,vp,valg,valpha,(int)(vx0),vis_integer,vdist,vrd_p);
									sample.AverageGenerator();
									results.push_back(sample.outputmain());
								}
								else if(constcoef==0){
									Samples<int> sample(vN,vTv[ii],vL,vcv,vhv,vpv,valg,valpha,(int)(vx0),vis_integer,vdist,vrd_p);
									sample.AverageGenerator();
									results.push_back(sample.outputmain());
								}
							}
						}//end of getting important results
						cout<<'\r'<<"Outputing data now...."<<setw(30)<<' '<<endl;
						fout<<"Period"<<endl;
						for(ii=0;ii<vTv_l;ii++)
							fout<<vTv[ii]<<endl;
						fout<<endl<<endl;
						fout<<"TotalCost"<<endl;
						for(ii=0;ii<vTv_l;ii++)
								fout<<results[ii][0]<<endl;
						fout<<endl<<endl;
						fout<<"TotalHolding Cost"<<endl;
						for(ii=0;ii<vTv_l;ii++)
								fout<<results[ii][1]<<endl;
						fout<<endl<<endl;
						fout<<"TotalBacklogging Cost"<<endl;
						for(ii=0;ii<vTv_l;ii++)
								fout<<results[ii][2]<<endl;
						fout<<endl<<endl;
						fout<<"TotalDemand"<<endl;
						for(ii=0;ii<vTv_l;ii++)
								fout<<results[ii][3]<<endl;
						fout<<endl<<endl;
						fout<<"TotalOrdering"<<endl;
						for(ii=0;ii<vTv_l;ii++)
								fout<<results[ii][4]<<endl;
						fout<<endl<<endl;
						fout<<"TotalCost/TotalDemand"<<endl;
						for(ii=0;ii<vTv_l;ii++)
								fout<<results[ii][0]/results[ii][0]<<endl;
						fout<<endl<<endl;
						

						fout.close();
						cout<<filename<<" data output finished."<<endl<<endl;
						break;
					}//end of open output file successfully
				}
			}//end of datatype 2
		}//end of for loop for different data sets
	}//end of open input
	return 1;
}




#endif