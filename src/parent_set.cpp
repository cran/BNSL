#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <vector>
#include <map>

#include "gc.h"
#include "parent_set.h"

inline int zerocount2(uint64_t x)
{
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >> 16);
    x = x | (x >> 32);
    return 64 - __builtin_popcountll( ~x );
}

// [[Rcpp::export]]
DataFrame parent(NumericMatrix df0, int h, int tw=0, int proc=0){
	int N=df0.ncol(), NN=N-1, n=df0.nrow(),i,j,k;
	if(tw==0)tw=NN;
	NumericMatrix df(n,N); for(j=0; j<h; j++)df(_,j)=df0(_,j); for(j=h; j<NN; j++)df(_,j)=df0(_,j+1); df(_,NN)=df0(_,h);
	Function co("combn"); 
	IntegerVector cho(tw+1);
	int LL=pow(NN,tw)+1;
	IntegerVector x(LL); IntegerMatrix y(LL,NN), w(LL,NN); NumericVector z(LL);
	int L=0, M=0;
	for(j=0; j<=tw; j++){
		IntegerMatrix u=co(NN,j);
		cho(j)=u.ncol();
		cho(0)=1;
		for(k=0; k<cho(j);k++){
			x(L+k)=0; z(L+k)=-1000000;
			IntegerVector uu=u(_,k);
			if(L+k>0){
				IntegerMatrix v=co(NN,j-1);
				IntegerMatrix uuu(j-1,j);
				for(i=0; i<j; i++){
					for(int r=0; r<i; r++)uuu(r,j-i-1)=uu(r);
					for(int r=i+1; r<j; r++)uuu(r-1,j-i-1)=uu(r);
				}
				i=0;
				for(int r=0; r<j; r++)
				{
					IntegerVector a=uuu(_,r), b=v(_,i); while(!setequal(a,b)){i++; b=v(_,i);}
					if(x(M+i)==1)x(L+k)=1; 
					if(z(L+k)<z(M+i)){y(L+k,_)=y(M+i,_); z(L+k)=z(M+i);};
					i++;					
				}
			}
			for(i=0; i<NN; i++)w(L+k,i)=0; for(i=0; i<j; i++)w(L+k,uu(i)-1)=1;
			if(x(L+k)==0){
				NumericMatrix df2(n,N);
				int r=0; for(i=0; i<NN; i++)if(w(L+k,i)==1)df2(_,r++)=df(_,i); df2(_,r)=df(_,NN);
				int m=table(df(_,NN)).size();
				IntegerMatrix T=fftable(df2(_,Range(0,r)), m); 
				if(z(L+k)>bound(T,m,proc))x(L+k)=1;
				else{double s=Bayes_score(T,m,proc); if(s>z(L+k)){y(L+k,_)=w(L+k,_); z(L+k)=s;}}
			}
		}
		M=L; L=L+cho(j);
	}
	DataFrame results = DataFrame::create(Named("w")=w,Named("y")=y,Named("z")=z); return(results);
}

// [[Rcpp::export]]
IntegerMatrix fftable(NumericMatrix df, int w){
    int u = df.nrow(), v = df.ncol();
    IntegerMatrix T(u,w);
    std::map<uint64_t, std::map<int, int> > counter;

    IntegerMatrix M(u, v - 1);
    IntegerVector S(v);
    S[0] = 0;

    for (int j = 0; j < v - 1; ++j) {
        std::map<int, int> val_map;
        int c = 0;
        for (int i = 0; i < u; ++i) {
            if (val_map.find(df(i, j)) == val_map.end()) { // not found
                val_map[df(i, j)] = c;
                ++c;
            }
            M(i, j) = val_map[df(i, j)];
        }
        S[j + 1] = S[j] + zerocount2(c - 1);
        if (S[j + 1] > 64) {
            stop("S cannot exceed 64!\n");
        }
    }

    for (int i = 0; i < u; ++i) {
        uint64_t b = 0;
        for (int j = 0; j < v - 1; ++j) {
            b |= (static_cast<uint64_t>(M(i, j)) << S[j]);
        }
        int val = df(i, v - 1);
        counter[b][val] += 1;
    }
    int row = 0;
    std::map<uint64_t, std::map<int, int> >::iterator itor1 = counter.begin();
    while (itor1 != counter.end()) {
        int col = 0;
        std::map<int, int>::iterator itor2 = itor1->second.begin();
        while (itor2 != itor1->second.end()) {
            T(row, col) = itor2->second;
            ++col;
            ++itor2;
        }
        ++row;
        ++itor1;
    }
    return T(Range(0, row - 1), _);
}


// [[Rcpp::export]]
double Bayes_score(IntegerMatrix T, int m, int proc=1){
	if(proc==1)return(Jeffreys_score(T,m));
	else if(proc==2)return(MDL_score(T,m));
	else if(proc==3)return(BDeu_score(T,m));
	else return(Jeffreys_score(T,m));	
}

// [[Rcpp::export]]
double bound(IntegerMatrix T, int m, int proc=1){
	if(proc==1)return(Jeffreys_bound(T,m));
	else if(proc==2)return(MDL_bound(T,m));
	else if(proc==3)return(BDeu_bound(T,m));
	else return(Jeffreys_bound(T,m));	
}


// [[Rcpp::export]]
double Jeffreys_score(IntegerMatrix T, int m){
        double s=0;
        int j, w=T.nrow();
        for(j=0; j<w; j++)s=s-gc(sum(T(j,_)),m*0.5)+gc_all(T(j,_),0.5);
        return (s);
}

// [[Rcpp::export]]
double BDeu_score(IntegerMatrix T, int m){
        double s=0;
        int j, w=T.nrow();
        for(j=0; j<w; j++)s=s-gc(sum(T(j,_)),1./w)+gc_all(T(j,_),1./m/w);
        return (s);
}

// [[Rcpp::export]]
double MDL_score(IntegerMatrix T, int m){
        double s=0;
        int j, w=T.nrow();
		int n=0; for(j=0; j<w; j++)n=n+sum(T(j,_));
		for(j=0; j<w; j++){double n_s=sum(T(j,_)); for(int k=0; k<m; k++)s=s+T(j,k)*log(T(j,k)/n_s); s=s-0.5*(m-1)*log(1.0*n);} 
        return (s);
}

// [[Rcpp::export]]
double Jeffreys_bound(IntegerMatrix T, int m){
 	double s=0;
 	int j, w=T.nrow();
 	for(j=0; j<w; j++)s=s+gc(sum(T(j,_)),0.5)-gc(sum(T(j,_)),0.5*m);
 	return (s);
}

// [[Rcpp::export]]
double MDL_bound(IntegerMatrix T, int m){
        int j, w=T.nrow();
		int n=0; for(j=0; j<w; j++)n=n+sum(T(j,_));
		double s=-0.5*(m-1)*w*log(1.0*n);
        return (s);
}

// [[Rcpp::export]]
double BDeu_bound(IntegerMatrix T, int m){
        int w=T.nrow();
		double s=pow(m,w);
        return (s);
}



