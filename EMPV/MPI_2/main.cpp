#define USE_MATH_DEFINES
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <functional>
#include <mpi.h>
#include <map>

using namespace std;

enum TTargetFunctionType {
	TFT_NotDefined = -1,
	TFT_Square,
	TFT_Rozenbrooke,
	TFT_Rastrigin
};

enum TCrossoverFunctionType {
	CFT_NotDefined = -1,
	CFT_Default,
	CFT_TwoPoint,
	CFT_OnePoint
};



class MyRandom {
public:

	MyRandom(){
		en = std::mt19937(rd());
		dist = uniform_real_distribution<double>(-100.0, 100.0);
	}

	double Generate(){
		return dist(en);
	}
    
 private:
 	std::random_device rd;
 	std::mt19937 en;
 	std::uniform_real_distribution<double> dist;
};

double frand()
{
	static MyRandom rnd;
	return rnd.Generate();
}

void shuffle( vector< vector< double> >& a )
{
	for( int k = 0; k < a.size(); k++ ){
		int l = rand() % a.size();
		for( int i = 0; i < a[k].size(); i++ ){
			swap( a[k][i], a[l][i] );
		}
	}
}

std::function< double(const vector<double>&)> getTargetFunc( TTargetFunctionType tft ){
	switch( tft ){
		case TFT_Square:
			return []( const vector<double>& a ){ 
				double sum = 0;
				for( double x : a ) {
					sum += x * x;
				}
				return sum;
			};
		case TFT_Rastrigin:
			return []( const vector<double>& a ){ 
				double sum = 10 * a.size();
				for( double x : a ) {
					sum += x * x - 10 * cos( 2 * M_PI * x );
				}
				return sum;
			};
		case TFT_Rozenbrooke:
			return []( const vector<double>& a ){ 
				double sum = 0;
				for( int i = 0; i < int( a.size() ) - 1; i++ ){
					double val  = ( a[i] * a[i] - a[i + 1]);
					sum += 100.0 * val * val + ( a[i] - 1.0 ) * ( a[i] - 1.0 );
				}
				return sum;
			};
		default:
			return []( const vector<double>& a ){ return 0; };
	}
}

std::function< void( vector< vector<double> >& )> getCrossoverFunc( TCrossoverFunctionType cft ){
	switch(cft){
		case CFT_Default:
			return []( vector< vector<double> >& P ){ 
						shuffle( P );
						for( int k = 0; k < P.size() / 2; k++ )
						{
							int a = 2 * k;
							int b = 2 * k + 1;
							
							for( int i = 0; i < P[a].size(); i++ ) {
								bool rnd = frand() < 0.0;
								if( rnd )
									swap( P[a][i], P[b][i] );
							}
						}
					};
		case CFT_OnePoint:
			return []( vector< vector<double> >& P ){ 
						shuffle( P );
						for( int k = 0; k < P.size() / 2; k++ )
						{
							int a = 2 * k;
							int b = 2 * k + 1;
							int j = rand() % P[a].size();

							for( int i = j; i < P[a].size(); i++ )
								swap( P[a][i], P[b][i] );
						}
					};
		case CFT_TwoPoint:
			return []( vector< vector<double> >& P ){ 
						shuffle( P );
						for( int k = 0; k < P.size() / 2; k++ )
						{
							int a = 2 * k;
							int b = 2 * k + 1;
							int j = rand() % P[a].size();
							int jj = rand() % P[a].size();
							if( jj < j ){
								swap(jj, j);
							}
							for( int i = j; i < jj; i++ )
								swap( P[a][i], P[b][i] );
						}
					};
	}
}






class Solution{
public:
	Solution(int n, int m, int T, TTargetFunctionType tft, int period, TCrossoverFunctionType cft, double _percentage) : 
		cnt(0), 
		P( vector< vector< double > >(n, vector< double > (m) ) ),
		recvbuf( vector< vector< double > >( static_cast<int>( _percentage * double(n) ) , vector< double >( m ) ) ),
		targetFunc( getTargetFunc( tft ) ),
		crossoverFunc( getCrossoverFunc( cft ) )
	{

		int rank, size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		const int totalSpecies = size * n;

		// cout << recvbuf.size() << endl;

		init( P );
		for( int t = 0, cnt = 0; t < T; cnt++, t++ )
		{
			select( P );
			crossoverFunc( P );
			mutate( P );
			if( cnt == period ){
				cnt = 0;
				migrate( P );
			}
			// printthebest(P)
			{
				double err = 0;
				double best = 1e18;
				for( const auto& v : P ){
					double f = targetFunc( v );
					double af = fabs( f );
					err += af;
					best = min( best, af );
				}
				double globalErr = 0;
				double globalBest = 0;
				MPI_Reduce( &err, &globalErr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Reduce( &best, &globalBest, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
				if ( rank == 0 ){
					globalErr /= totalSpecies;
					cout << "iter " << t << " avgErr " << globalErr << " bestErr " << globalBest << endl;	
				}
			}
		}

		// printthebest( P );
	}

	void init( vector< vector< double > >& a )
	{
		for( int k = 0; k < a.size(); k++ )
			for( int i = 0; i < a[k].size(); i++ )
				a[k][i] = frand();
	}

	void select( vector< vector< double > >& P )
	{
		double pwin = 0.75;
		shuffle( P );
		for( int k = 0; k < P.size() / 2; k++ )
		{
			int a = 2 * k;
			int b = 2 * k+1;
			double fa = targetFunc( P[a] );
			double fb = targetFunc( P[b] );
			double p = ( frand() + 100 ) / 200;
			if( fa < fb && p < pwin || fa > fb && p > pwin ){
				P[b] = P[a];
			}
			else {
				P[a] = P[b];
			}
		}
	}
	void mutate( vector< vector< double > >& P )
	{
		double pmut = 0.1;
		static double mult = 0.1;
		mult = mult * 0.9999;
		for( int k = 0; k < P.size(); k++ )
			if( ( frand() + 100 ) / 200 < pmut )
				for( int i = 0; i < P[k].size(); i++ )
					P[k][i] += frand() * mult;
	}

	void printthebest( vector< vector< double > >& P )
	{
		double k0 = -1;
		double f0 = -1;
		for( int k = 0; k < P.size(); k++ )
		{
			double f = targetFunc( P[k] );
			if( f < f0 || k == 0 )
			{
				f0 = f;
				k0 = k;
			}
		}
		cout << f0 << ": ";
		for( int i = 0; i < P[k0].size(); i++ )
			cout << P[k0][i] << " ";
		cout << endl;
	}

	void migrate( vector< vector<double> >& sendbuf ){
		
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		
		const int sendTo = (rank + 1) % size;
		const int recvFrom = (size + rank - 1) % size;

		const int part = recvbuf.size();
		vector< MPI_Request > rq( part * 2 );
		// if( rank == 0 ) {
		// 	// cout << "part " << part << endl;
		// 	for( int i =0 ; i < part; i++ ){
		// 		cout << sendbuf[i].size() << " " << sendbuf[i].data() << endl; 
		// 	}
		// }

		for( int i = 0; i < part; i++ ){
			MPI_Isend( sendbuf[i].data(), sendbuf[i].size(), MPI_DOUBLE, sendTo, i, MPI_COMM_WORLD, &rq[i]);
		}

		for( int i = 0; i < part; i++ ){
			MPI_Irecv( recvbuf[i].data(), recvbuf[i].size(), MPI_DOUBLE, recvFrom, i, MPI_COMM_WORLD, &rq[part + i]);
		}
		
		MPI_Waitall( 2 * part, rq.data(), MPI_STATUSES_IGNORE);
		
		for( int i = 0; i < part; i++ ){
			sendbuf[i] = recvbuf[i];
		}
		//cout << "here";
	}

private:
	int cnt;
	vector< vector< double > > P;
	vector< vector< double > > recvbuf;
	function< double(const vector<double>& ) > targetFunc;
	function< void( vector< vector<double> >& ) > crossoverFunc;
};



int main( int argc, char** argv )
{
	MPI_Init( &argc, &argv );
	if( argc < 8 ){
		cout << "use ./main n m T crossoverTypeString targetFunctionTypeString period migrationPercentage " << endl;
		cout << "crossoverTypeString: default, twopoint, onepoint " << endl;
		cout << "targetFunctionTypeString: square, rastrigin, rozenbrooke " << endl;
		return 0;
	}


	const int n = atoi( argv[1] ); // 10
	const int m = atoi( argv[2] ); // 20
	const int T = atoi( argv[3] ); // 10
	const string crossoverType( argv[4] );
	const string targetFunctionTypeString( argv[5] );
	const int period = atoi( argv[6] );
	const double perc = atof( argv[7] );
	map< string, TTargetFunctionType > targetFuncType;
	targetFuncType[ "square" ] = TFT_Square;
	targetFuncType[ "rastrigin" ] = TFT_Rastrigin;
	targetFuncType[ "rozenbrooke" ] = TFT_Rozenbrooke;

	map< string, TCrossoverFunctionType > crossFuncType;
	crossFuncType[ "default" ] = CFT_Default;
	crossFuncType[ "onepoint" ] = CFT_OnePoint;
	crossFuncType[ "twopoint" ] = CFT_TwoPoint;
	
	Solution( n, m, T, targetFuncType[targetFunctionTypeString],
			 period, crossFuncType[crossoverType], perc );

	MPI_Finalize();
	return 0;
} 
