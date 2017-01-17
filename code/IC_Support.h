#ifndef IC_SUPPORT_H
#define IC_SUPPORT_H
#include <cmath>
#include <time.h>
using namespace std;
const double PI=3.14159265359;


//standard distribution sample function
double Unif();
int Bernoulli(double p);
int Binomial(int n,double p);
double Normal(double mean,double std);
double Gamma();

double Unif(){
	//srand(time(NULL));
	double x = ((double) rand() / (RAND_MAX));
	return x;
}

int Bernoulli(double p){
	double x=Unif();
	if(x>=p) return 0;
	else return 1;
}

int Binomial(int n,double p){
	int temp=0,i;
	for(i=1;i<=n;i++)
		temp=temp+Bernoulli(p);
	return temp;
}

double Normal(double mean,double std){
	double u,v;
	u=Unif();
	v=Unif();
	return sqrt(-2 * log(u)) * cos(2 * PI * v) * std + mean;;

}








void Explanation(){
	cout<<"Welcome to use the stochastic inventory control simulation program. The copyright of this program belongs to Blake, Hao and Feng."<<endl;
}

void DC_Exp(){
	cout<<"Welcome to use the DataCollector of the stochastic inventory control."<<endl;
	cout<<"In your input file, there should be just one number in the first line specify how many data sets you are going to generate. ";
	cout<<"Every data set should be separated by an empty new line.";
	cout<<endl;
	cout<<"For every data set, please specify the following parameters: "<<endl;
	cout<<" 0. You want to collect averaging data of sample path for given parameters OR you want to collect data for a range of different PLANNING HORIZON. ";
	cout<<"(Put 1 if you want to do the first, 2 for the second.)"<<endl;
	cout<<endl;

	cout<<"If your answer is 1 for parameter 0, then specify the followings:"<<endl;
	cout<<" 1. How many iterations you want to compute to get the average data;(int)"<<endl;
	cout<<" 2. What is your planning horizon;(int)"<<endl;
	cout<<" 3. What is your leading time;(int)"<<endl;
	cout<<" 4. What is your discount rate and initial inventory;(alpha double, between 0 and 1, inventory double)"<<endl;
	cout<<" 5. Is this a integer typed problem or not;(1 for yes, 0 for no)"<<endl;
	cout<<" 6. What is the algorithm; (1: DP, 2: Myopic, 3: Dual Balancing)"<<endl;
	cout<<" 7. What is your distribution type? (1: Binomial, 2: AR1 Bernoulli, 3: ind Normal, 4: AR1 Normal, 5: ARn Normal)"<<endl;
	cout<<" 8. What are your distribution parameters?"<<endl;
	cout<<" 9. Whether your cost coefficient are going to be constant; (1 for constant, 0 for vector)"<<endl;
	cout<<" 10,11,12. Specify three vectors or three numbers depending on 9.(order is buying cost, holding cost then backlogging cost.)"<<endl;
	cout<<endl;
	cout<<"If your answer is 2 for parameter 0, then for parameter 2(planning horizon), specify a number indicate how many different horizons, and follow by a vector"<<endl;
	cout<<endl;

}

#endif