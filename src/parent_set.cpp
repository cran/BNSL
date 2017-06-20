#include <Rcpp.h>
using namespace Rcpp;
DataFrame parent(NumericMatrix, int, int, int);
IntegerMatrix fftable(NumericMatrix, int);
double gc(int, double);
double gc_all(IntegerVector, double);
double Bayes_score(IntegerMatrix, int, int, double, int, int);
double quotient_Jeffreys_score(IntegerMatrix, int, double, int);
double Jeffreys_score(IntegerMatrix, int);
double MDL_score(IntegerMatrix, int, double, int);
double BDeu_score(IntegerMatrix, int, int);
double bound(IntegerMatrix, int, int, int, int);
double quotient_Jeffreys_bound(IntegerMatrix, int, int, int);
double Jeffreys_bound(IntegerMatrix, int);
double MDL_bound(IntegerMatrix, int, int, int);
double BDeu_bound(IntegerMatrix, int);

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
        int N=df0.ncol(), NN=N-1, n=df0.nrow(), L=pow(2,NN),i,j,k;
		if(tw==0)tw=NN;
        NumericMatrix df(n,N);
        for(j=0; j<h; j++)df(_,j)=df0(_,j); for(j=h; j<NN; j++)df(_,j)=df0(_,j+1); df(_,NN)=df0(_,h);
        IntegerVector x(L), y(L); NumericVector z(L);
		for(i=0; i<L; i++)z(i)=-1000000;
        for(i=0; i<L; i++){
                IntegerVector set(NN);
				j=i; for(k=0; k<NN; k++)if(j>=pow(2,NN-1-k)){j=j-pow(2,NN-1-k); set(NN-1-k)=1;} else set(NN-1-k)=0;
                for(k=0; k<NN; k++)if(set(k)==1){
                        j=i-pow(2,k);
						if(x(j)==1)x(i)=1;
                        if(z(i)<z(j)){y(i)=y(j); z(i)=z(j);}
                }
                if(sum(set)<=tw&&x(i)==0){
                        NumericMatrix df2(n,N);
                        j=0; for(k=0; k<NN; k++)if(set(k)==1)df2(_,j++)=df(_,k); df2(_,j)=df(_,NN);
						int ss=1; for(k=0; k<j; k++)ss=ss*table(df2(_,k)).size();
                        int m=table(df(_,NN)).size();
						IntegerMatrix T=fftable(df2(_,Range(0,j)), m);
						double bb=bound(T,m,proc,n,ss);
                        if(z(i)>bb)x(i)=1;
                        else{double s=Bayes_score(T,m,proc,bb,n,ss); if(s>z(i)){y(i)=i; z(i)=s;}}
                }
        }
		DataFrame results = DataFrame::create(Named("y")=y,Named("z")=z); return(results);
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
double Bayes_score(IntegerMatrix T, int m, int proc=0, double s=0, int n=0, int ss=1){
	if(proc==0)return(quotient_Jeffreys_score(T,m,s,n));
	else if(proc==1)return(Jeffreys_score(T,m));	
	else if(proc==2)return(MDL_score(T,m,s,n));
	else if(proc==3)return(BDeu_score(T,m,ss));
	else return(quotient_Jeffreys_score(T,m,s,n));	
}

double quotient_Jeffreys_score(IntegerMatrix T, int m, double s, int n){
        int w=T.nrow();
        for(int j=0; j<w; j++)s=s-gc(sum(T(j,_)),0.5)+gc_all(T(j,_),0.5);
        return (s);
}

// [[Rcpp::export]]
double Jeffreys_score(IntegerMatrix T, int m){
        int w=T.nrow();
        double s=0;
        for(int j=0; j<w; j++)s=s-gc(sum(T(j,_)),m*0.5)+gc_all(T(j,_),0.5);
        return (s);
}

// [[Rcpp::export]]
double MDL_score(IntegerMatrix T, int m, double s, int n){
        int w=T.nrow();
		for(int j=0; j<w; j++){
			double n_s=sum(T(j,_)); 
			for(int k=0; k<m; k++)s=s+T(j,k)*log(T(j,k)/n_s);
		} 
        return (s);
}

// [[Rcpp::export]]
double BDeu_score(IntegerMatrix T, int m, int ss){
        int w=T.nrow();
        double s=0;
        for(int j=0; j<w; j++)s=s-gc(sum(T(j,_)),1./ss)+gc_all(T(j,_),1./m/ss);
        return (s);
}

// [[Rcpp::export]]
double bound(IntegerMatrix T, int m, int proc=0, int n=0, int ss=1){
	if(proc==0)return(quotient_Jeffreys_bound(T,m,n,ss));
	else if(proc==1)return(Jeffreys_bound(T,m));
	else if(proc==2)return(MDL_bound(T,m,n,ss));
	else if(proc==3)return(BDeu_bound(T,m));
	else return(Jeffreys_bound(T,m));
}

// [[Rcpp::export]]
double Jeffreys_bound(IntegerMatrix T, int m){
		int w=T.nrow();
		double s=0;
		for(int j=0; j<w; j++)s=s+gc_all(T(j,_),0.5)-gc_all(T(j,_),0.5*m);
		return (s);
}

// [[Rcpp::export]]
double quotient_Jeffreys_bound(IntegerMatrix T, int m, int n, int ss){
		return(gc(n,0.5*ss)-gc(n,0.5*ss*m));
}

// [[Rcpp::export]]
double MDL_bound(IntegerMatrix T, int m, int n, int ss){
        return (-0.5*(m-1)*ss*log(1.0*n));
}

// [[Rcpp::export]]
double BDeu_bound(IntegerMatrix T, int m){
		int w=T.nrow();
        return (-w*log(1.0*m));
}



