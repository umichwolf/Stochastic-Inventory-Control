#ifndef IC_RANDOM_H
#define IC_RANDOM_H
#include"IC_Policy.h"
#include"IC_Support.h"
using namespace std;

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


template<class W>
vector<double> ConditionalExp(double q,Policy<W> &policy){
	// cout << "Here!" << endl;
	//cout<<"Alert!!!!!!!!!!!!!!!! The ConditionalExp is called!!!!!!!!!!"<<endl;
	if(policy.dist==1)
		return BinomialCE(q,policy);
	else if(policy.dist==2)
		return IncreBinomialCE(q,policy);
	else if(policy.dist==3)
		return NormalCE(q,policy);
	else if (policy.dist == 4)
		return RA1NormalCE(q,policy);
	else if (policy.dist == 5)
		return RAnNormalCE(q,policy);
	vector<double> temp(3,1);
	return temp;
}

struct X // helper functions for compute statistics about normal dist
{
	static double phi(double x) { //Pdf of normal dist
		return 1/sqrt(2*PI) * exp(-0.5 * x*x); 
	}

	static double Phi(double x) // A numatrical method to compute cdf of standard normal distribution
	{
	    // constants
	    double a1 =  0.254829592;
	    double a2 = -0.284496736;
	    double a3 =  1.421413741;
	    double a4 = -1.453152027;
	    double a5 =  1.061405429;
	    double p  =  0.3275911;
	 
	    // Save the sign of x
	    int sign = 1;
	    if (x < 0)
	        sign = -1;
	    x = fabs(x)/sqrt(2.0);
	 
	    // A&S formula 7.1.26
	    double t = 1.0/(1.0 + p*x);
	    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
	 
	    return 0.5*(1.0 + sign*y);
	}

	static double exp_tranc_normal(double mu, double sigma, double a, double b) { 
		double alpha = (a - mu) / sigma;
		double beta = (b - mu) / sigma;
		double Z = Phi(beta) - Phi(alpha);
		//cout<<"Z:"<<Z<<' '<<phi(alpha) - phi(beta)<<" return:"<<temp<<endl;
		return mu*Z + (phi(alpha) - phi(beta))*sigma;
	}


	static double c_mean(double d, double inc,double mu, double sigma) { // compute mean of RA1 under normal pertbation
		return (inc+1)*d + (inc + 1)*(inc + 2) / 2.0 * mu;
	}

	static double c_sigma(double d, double inc, double mu, double sigma) { // compute sigma of RA1 normal
		double sum = 0;
		for (int i = 1; i<=inc+1; i++) {
			sum += i * i * sigma * sigma;
		}
		return sqrt(sum);
	}

	static vector< vector<double> > mean_std(int n,vector<double> d, double inc, double mu, double sigma) {
		vector< vector<double> > val;
		vector< vector<double> > sum;
		vector<double> temp(n+inc+2,0);
		int i,j,k;
		for(i=0;i<=inc+1+n;i++){
			val.push_back(temp);
			sum.push_back(temp);
		}

		for(i=1;i<=n+1+inc;i++)
			val[i][i]=1;

		for(i=n+1;i<=inc+1+n;i++){
			for(j=1;j<=i-1;j++){
				for(k=1;k<=n;k++){
					val[i][j]=val[i][j]+val[i-k][j]/n;
				}
			}
		}
		//cout<<n<<' '<<d.size()<<' '<<inc<<' '<<mu<<' '<<sigma<<endl;

		for(i=n+1;i<=inc+1+n;i++){
			for(j=1;j<=inc+1+n;j++){
				sum[i][j]=sum[i-1][j]+val[i][j];
			}
		}
		
		vector< vector<double> > ret;
		vector<double> tempret;
		double mean,stdd;

		for(i=n+1;i<=inc+1+n;i++){
			tempret.clear();
			mean=0,stdd=0;
			for(j=1;j<=n;j++)
				mean=mean+sum[i][j]*d[j];
			for(j=n+1;j<=i;j++){
				mean=mean+sum[i][j]*mu;
				stdd=stdd+sum[i][j]*sum[i][j]*sigma*sigma;
			}
			tempret.push_back(mean);
			tempret.push_back(sqrt(stdd));
			ret.push_back(tempret);
		}
		return ret;
	}


	// static vector<double> (double p, double n) { // Probability masses for Bernulli
	// 	vector<double> v (n+1,0);
	// 	for (int i = 0; i <= n; i++) {
	// 		v[i] = 
	// 	}
	// }
};

template<class W>
int RandomDemand(Policy<W> &policy){
	//cout<<"Alert!!!!!!!!!!!!!!!! The RandomDemand is called!!!!!!!!!!"<<endl;
	if(policy.dist==1)
		policy.d[policy.cur]=Binomial((int)policy.rd_p[0],policy.rd_p[1]);
	else if(policy.dist==2) {
		if (policy.cur == 1) 
			policy.d[policy.cur] = policy.rd_p[2] + Binomial(policy.rd_p[0],policy.rd_p[1]);
		else 
			policy.d[policy.cur]= policy.d[policy.cur - 1] + Binomial(policy.rd_p[0],policy.rd_p[1]);
	}
		// policy.d[policy.cur]=policy.d[policy.cur]+Binomial((int)policy.rd_p[0],policy.rd_p[1]);
	else if(policy.dist==3)
		policy.d[policy.cur]=Normal(policy.rd_p[0],policy.rd_p[1]);
	else if (policy.dist==4) {
			if (policy.cur == 1) 
				policy.d[policy.cur] = policy.rd_p[2] + Normal(policy.rd_p[0],policy.rd_p[1]);
			else 
				policy.d[policy.cur]= policy.d[policy.cur - 1] + Normal(policy.rd_p[0],policy.rd_p[1]);
		}
	else if(policy.dist==5){
		double d=0;
		double n=policy.rd_p[3];
		if(policy.cur-n<1){
			if(policy.cur == 1)
				d=policy.rd_p[2];
			else{
				for(int i=1;i<=policy.cur-1;i++)
					d=d+policy.d[i];
				d=d/(policy.cur-1);
			}
		}
		else{
			for(int i=policy.cur-n;i<=policy.cur-1;i++)
				d=d+policy.d[i];
			d=d/n;
		}
		policy.d[policy.cur]=d+Normal(policy.rd_p[0],policy.rd_p[1]);
	}
	return 1;
}

template<class W> 
vector<double> BinomialCE(double q,Policy<W> &P){
	struct Y {
			static vector<double> DualBalancing(double q, Policy<W> &P) {
			vector<double> v;
			double mu = P.rd_p[0] * P.rd_p[1];
			double sigma = sqrt(P.rd_p[0] * (1 - P.rd_p[1])*P.rd_p[1]);
			double sum = 0;
			double h_hat;
			for (int j = P.cur+P.L; j <= P.T; j++) {
				if(j == P.T) {
					h_hat = P.h[j] + P.c[j];
				} else {
					h_hat = P.h[j] + P.c[j] - P.c[j+1];
				}
				double inc = h_hat * X::exp_tranc_normal(q - ((j-P.cur+1)*mu - P.x[P.cur-1]), sqrt(j - P.cur + 1.0)*sigma,0,q);
				if (inc == 0) break;
				sum += inc;
			}
			v.push_back(sum);
			double p_bar;
			if (P.cur == P.T)
				p_bar = P.p[P.cur] - P.c[P.cur];
			else
				p_bar = P.p[P.cur] - P.c[P.cur] + P.c[P.cur+1];
			v.push_back(p_bar*X::exp_tranc_normal((P.L+1)*mu - P.x[P.cur-1] - q,sqrt(P.L+1.0)*sigma,0,100000000000));
			return v;
		}
		static vector<double> Myopic(double q, Policy<W> &P) {
			vector<double> v;
			double mu = P.rd_p[0] * P.rd_p[1];
			double sigma = sqrt(P.rd_p[0] * (1 - P.rd_p[1])*P.rd_p[1]);
			double h_hat, p_hat;
			if (P.cur + P.L == P.T) {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L];
			} else {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L] - P.c[P.cur + P.L + 1];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L] + P.c[P.cur + P.L + 1];
			}
			double result = h_hat * X::exp_tranc_normal(P.x[P.cur-1] + q - (P.L+1)*mu, sqrt(P.L+1.0)*sigma ,0,100000000000) +
						   p_hat * X::exp_tranc_normal(-P.x[P.cur-1] - q + (P.L+1)*mu, sqrt(P.L+1.0)*sigma ,0,100000000000);
			v.push_back(result);
			return v;
		}
	};
	vector<double> v;
	if(P.alg == 2) v  = Y::Myopic(q,P);
	else if(P.alg == 3) v = Y::DualBalancing(q,P);
	return v;
}

template<class W> 
vector<double> IncreBinomialCE(double q,Policy<W> &P){
	struct Y {
		static vector<double> DualBalancing(double q, Policy<W> &P) {
			double mu = P.rd_p[0] * P.rd_p[1];
			double sigma = sqrt(P.rd_p[0] * (1 - P.rd_p[1])*P.rd_p[1]);
			double L = P.L;
		    double d = P.d[P.cur - 1];
		    if(P.cur == 1) d = P.rd_p[2];
		    double x = P.x[P.cur - 1];
		    double y = x + q;
		    double h_hat;
			vector<double> v;
			double sum = 0;
			for (int j = P.cur+P.L; j <= P.T; j++) {
				if(j == P.T) {
					h_hat = P.h[j] + P.c[j];
				} else {
					h_hat = P.h[j] + P.c[j] - P.c[j+1];
				}
				double inc = h_hat * X::exp_tranc_normal(q - (X::c_mean(d,j-P.cur,mu,sigma) - x), X::c_sigma(d,j-P.cur,mu,sigma),0,q);
				if (inc == 0) break;
				sum += inc;
			}
			v.push_back(sum);
			double p_bar;
			if (P.cur == P.T)
				p_bar = P.p[P.cur] - P.c[P.cur];
			else
				p_bar = P.p[P.cur] - P.c[P.cur] + P.c[P.cur+1];
			v.push_back(p_bar*X::exp_tranc_normal(X::c_mean(d,L,mu,sigma) - y,X::c_sigma(d,L,mu,sigma),0,100000000000));
			return v;
		}
		static vector<double> Myopic(double q, Policy<W> &P) {
			double mu = P.rd_p[0] * P.rd_p[1];
			double sigma = sqrt(P.rd_p[0] * (1 - P.rd_p[1])*P.rd_p[1]);
			double L = P.L;
		    double d = P.d[P.cur - 1];
		    if(P.cur == 1) d = P.rd_p[2];
		    double x = P.x[P.cur - 1];
		    double y = x + q;
			vector<double> v;
			double h_hat, p_hat;
			if (P.cur + P.L == P.T) {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L];
			} else {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L] - P.c[P.cur + P.L + 1];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L] + P.c[P.cur + P.L + 1];
			}
			double result = h_hat*X::exp_tranc_normal(y - X::c_mean(d,L,mu,sigma),X::c_sigma(d,L,mu,sigma),0,100000000000) +
						    p_hat*X::exp_tranc_normal(X::c_mean(d,L,mu,sigma) - y,X::c_sigma(d,L,mu,sigma),0,100000000000);
			v.push_back(result);
			return v;
		}
	};

	vector<double> v;
	if(P.alg == 2) v  = Y::Myopic(q,P);
	else if(P.alg == 3) v = Y::DualBalancing(q,P);
	return v;
}



template<class W>
vector<double> NormalCE(double q, Policy<W> &P){

	struct Y {
			static vector<double> DualBalancing(double q, Policy<W> &P) {
			vector<double> v;
			double mu = P.rd_p[0];
			double sigma = P.rd_p[1];
			double sum = 0;
			double h_hat;
			for (int j = P.cur+P.L; j <= P.T; j++) {
				if(j == P.T) {
					h_hat = P.h[j] + P.c[j];
				} else {
					h_hat = P.h[j] + P.c[j] - P.c[j+1];
				}
				double inc = h_hat * X::exp_tranc_normal(q - ((j-P.cur+1)*mu - P.x[P.cur-1]), sqrt(j - P.cur + 1.0)*sigma,0,q);
				if (inc == 0) break;
				sum += inc;
			}
			v.push_back(sum);
			double p_bar;
			if (P.cur == P.T)
				p_bar = P.p[P.cur] - P.c[P.cur];
			else
				p_bar = P.p[P.cur] - P.c[P.cur] + P.c[P.cur+1];
			v.push_back(p_bar*X::exp_tranc_normal((P.L+1)*mu - P.x[P.cur-1] - q,sqrt(P.L+1.0)*sigma,0,100000000000));
			return v;
		}
		static vector<double> Myopic(double q, Policy<W> &P) {
			vector<double> v;
			double mu = P.rd_p[0];
			double sigma = P.rd_p[1];
			double h_hat, p_hat;
			if (P.cur + P.L == P.T) {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L];
			} else {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L] - P.c[P.cur + P.L + 1];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L] + P.c[P.cur + P.L + 1];
			}
			double result = h_hat * X::exp_tranc_normal(P.x[P.cur-1] + q - (P.L+1)*mu, sqrt(P.L+1.0)*sigma ,0,100000000000) +
						   p_hat * X::exp_tranc_normal(-P.x[P.cur-1] - q + (P.L+1)*mu, sqrt(P.L+1.0)*sigma ,0,100000000000);
			v.push_back(result);
			return v;
		}
	};
	vector<double> v;
	if(P.alg == 2) v  = Y::Myopic(q,P);
	else if(P.alg == 3) v = Y::DualBalancing(q,P);
	return v;
}

template<class W> 
vector<double> RA1NormalCE(double q, Policy<W> &P){
	struct Y {
		static vector<double> DualBalancing(double q, Policy<W> &P) {
			double mu = P.rd_p[0];
			double sigma = P.rd_p[1];
			double L = P.L;
		    double d = P.d[P.cur - 1];
		    if(P.cur == 1) d = P.rd_p[2];
		    double x = P.x[P.cur - 1];
		    double y = x + q;
		    double h_hat;
			vector<double> v;
			double sum = 0;
			for (int j = P.cur+P.L; j <= P.T; j++) {
				if(j == P.T) {
					h_hat = P.h[j] + P.c[j];
				} else {
					h_hat = P.h[j] + P.c[j] - P.c[j+1];
				}
				double inc = h_hat * X::exp_tranc_normal(q - (X::c_mean(d,j-P.cur,mu,sigma) - x), X::c_sigma(d,j-P.cur,mu,sigma),0,q);
				if (inc == 0) break;
				sum += inc;
			}
			v.push_back(sum);
			double p_bar;
			if (P.cur == P.T)
				p_bar = P.p[P.cur] - P.c[P.cur];
			else
				p_bar = P.p[P.cur] - P.c[P.cur] + P.c[P.cur+1];
			v.push_back(p_bar*X::exp_tranc_normal(X::c_mean(d,L,mu,sigma) - y,X::c_sigma(d,L,mu,sigma),0,100000000000));
			return v;
		}
		static vector<double> Myopic(double q, Policy<W> &P) {
			double mu = P.rd_p[0];
			double sigma = P.rd_p[1];
			double L = P.L;
		    double d = P.d[P.cur - 1];
		    if(P.cur == 1) d = P.rd_p[2];
		    double x = P.x[P.cur - 1];
		    double y = x + q;
			vector<double> v;
			double h_hat, p_hat;
			if (P.cur + P.L == P.T) {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L];
			} else {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L] - P.c[P.cur + P.L + 1];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L] + P.c[P.cur + P.L + 1];
			}
			double result = h_hat*X::exp_tranc_normal(y - X::c_mean(d,L,mu,sigma),X::c_sigma(d,L,mu,sigma),0,100000000000) +
						    p_hat*X::exp_tranc_normal(X::c_mean(d,L,mu,sigma) - y,X::c_sigma(d,L,mu,sigma),0,100000000000);
			v.push_back(result);
			return v;
		}
	};

	vector<double> v;
	if(P.alg == 2) v  = Y::Myopic(q,P);
	else if(P.alg == 3) v = Y::DualBalancing(q,P);
	return v;
}

template<class W> 
vector<double> RAnNormalCE(double q, Policy<W> &P){
	struct Y {
		static vector<double> DualBalancing(double q, Policy<W> &P) {
			double mu = P.rd_p[0];
			double sigma = P.rd_p[1];
			double n=P.rd_p[3];
			double L = P.L;
		    vector<double> d(1,0);
		    double x = P.x[P.cur - 1];
		    double y = x + q;
		    double h_hat;
			vector<double> v;
			vector< vector<double> > ms;
			double sum = 0;

			if(P.cur-n<1){
				if(P.cur == 1)
					d.push_back(P.rd_p[2]);
				else{
					for(int i=1;i<=P.cur-1;i++)
						d.push_back(P.d[i]);
				}
				ms=X::mean_std(max(P.cur-1,1),d,P.T-P.cur,mu,sigma);
			}
			else{
				for(int i=P.cur-n;i<=P.cur-1;i++)
					d.push_back(P.d[i]);
				ms=X::mean_std(n,d,P.T-P.cur,mu,sigma);
			}
			for (int j = P.cur+P.L; j <= P.T; j++) {
				if(j == P.T) {
					h_hat = P.h[j] + P.c[j];
				} else {
					h_hat = P.h[j] + P.c[j] - P.c[j+1];
				}
				double inc = h_hat * X::exp_tranc_normal(q -( ms[j-P.cur][0] - x), ms[j-P.cur][1],0,q);
				// cout << "mean: " << ms[j-P.cur][0] << endl;
				// getchar();
				// getchar();
				if (inc == 0) break;
				sum += inc;
			}
			
			//cout<<"We are good here"<<P.cur<<endl;
			v.push_back(sum);
			double p_bar;
			if (P.cur == P.T)
				p_bar = P.p[P.cur] - P.c[P.cur];
			else
				p_bar = P.p[P.cur] - P.c[P.cur] + P.c[P.cur+1];
			v.push_back(p_bar*X::exp_tranc_normal(ms[P.L][0] - y,ms[P.L][1],0,100000000000));
			return v;
		}
		static vector<double> Myopic(double q, Policy<W> &P) {
			double mu = P.rd_p[0];
			double sigma = P.rd_p[1];
			double n=P.rd_p[3];
			double L = P.L;
		    vector<double> d(1,0);
		    double x = P.x[P.cur - 1];
		    double y = x + q;
		    double h_hat,p_hat;
			vector<double> v;
			vector< vector<double> > ms;
			double sum = 0;

			if (P.cur + P.L == P.T) {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L];
			} else {
				h_hat = P.h[P.cur + P.L] + P.c[P.cur+P.L] - P.c[P.cur + P.L + 1];
				p_hat = P.p[P.cur + P.L] - P.c[P.cur+P.L] + P.c[P.cur + P.L + 1];
			}

			if(P.cur-n<1){
				if(P.cur == 1)
					d.push_back(P.rd_p[2]);
				else{
					for(int i=1;i<=P.cur-1;i++)
						d.push_back(P.d[i]);
				}
				ms=X::mean_std(max(P.cur-1,1),d,P.L,mu,sigma);
			}
			else{
				for(int i=P.cur-n;i<=P.cur-1;i++)
					d.push_back(P.d[i]);
				ms=X::mean_std(n,d,P.L,mu,sigma);
			}


			double result = h_hat*X::exp_tranc_normal(y - ms[L][0],ms[L][1],0,100000000000) +
						    p_hat*X::exp_tranc_normal(ms[L][0] - y,ms[L][1],0,100000000000);
			v.push_back(result);
			return v;
		}
	};

	vector<double> v;
	if(P.alg == 2) v  = Y::Myopic(q,P);
	else if(P.alg == 3) v = Y::DualBalancing(q,P);
	return v;
}



#endif