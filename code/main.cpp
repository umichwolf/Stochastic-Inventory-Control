#include<vector>
#include<iostream>
//#include<cstdlib>
#include"IC_ParameterData.h"
#include"IC_Policy.h"
#include"IC_Random.h"
#include"IC_DualBalancing.h"
#include"IC_DP.h"
#include"IC_Myopic.h"
#include"IC_Statistics.h"
#include"IC_Support.h"
using namespace std;

int main(int argc,char **argv){
	srand(time(NULL));
	DataCollector();
	
	system("pause");
	return 0;
}

	/*
	int test[21]={0};
	int i,tempi;
	double temp;
	srand(time(NULL));
	for(i=0;i<=200;i++){
		temp=Normal(0,1);
		if(-10.0<=temp&&temp<10.0){
			tempi=(int)(temp+10);
		cout<<temp<<' '<<tempi<<' ';
			test[tempi]=test[tempi]+1;
		}
		else
			test[20]=test[20]+1;
	}
	cout<<endl<<endl;
	for(i=0;i<=20;i++)
		cout<<test[i]<<endl;
*/



/*

	vector<double> v(100,0.5);
	Policy<double> p1(100,2,v,v,v,3,0.99,10,0);
	RandomDemand(p1);
	vector<double> u;
	double y=10;
	u=ConditionalExp(y,p1);
	//cout<<u[0]<<endl<<u[2]<<endl<<u[3]<<endl;
	//p1.output();
	
	Samples<double> s1;
	Samples<int> s2(13,100,2,v,v,v,3,0.98,10,0);
	s2.AverageGenerator();
	s2.output();


*/




	
	/*
	ParameterData ic1;
	ic1.output();
	ic1.output_t(1);
	ParameterData ic2(5,2,0.5,0.3,0.2);
	ic2.output();

	Policy<double> db1;
	db1.output();
	db1.ParameterData::output();	
	
	
	vector<double> v;
	v.push_back(0.22);
	v.push_back(0.22);
	v.push_back(0.22);
	v.push_back(0.22);
	v.push_back(0.22);
	
	Policy<double> p1(5,2,v,v,v);
	p1.output();
	p1.ParameterData::output();
	*/