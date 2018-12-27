#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <bitset>
#include <chrono>
#include <deque>
#include <queue>
#include <fstream>
#include <stack>

using namespace std;

int INF = 1e9;

class Solution {
public:
	class Segment : public pair<int, int> {
	public:
		Segment() = default;
		Segment( int a, int b ) : pair<int, int>( a, b ) {}
		const int& L() const { return first; }
		int& L() { return first; }
		const int& R() const { return second; }
		int& R() { return second; }
		int Length() { return R() - L(); }
		const int Length() const { return R() - L(); }
	};

	static Segment Unite( const Segment& a, const Segment& b )
	{
		Segment res;
		res.L() = min( a.L(), b.L() );
		res.R() = max( a.R(), b.R() );
		if( res.Length() > a.Length() + b.Length() ) {
			res.L() = -1;
			res.R() = -1;
		}
		return res;
	}
	
	class SubtreeDescription : public vector<Segment> {
	public:
		size_t SegCnt() { return size(); }
		const size_t SegCnt() const { return size(); }
		int GetSubTreeSize()
		{
			if( !_upToDate ) {
				UpdateSubtreeSize();
			}
			return _subTreeSize;
		}
		void UpdateSubtreeSize()
		{
			_subTreeSize = 0;
			for( int i = 0; i < SegCnt(); i++ ) {
				_subTreeSize += ( *this )[i].Length();
			}
			_upToDate = true;
		}
		void SetOutdated()
		{
			_upToDate = false;
		}
	private:
		int _subTreeSize = 0;
		bool _upToDate = false;
	};

	static SubtreeDescription Unite( const SubtreeDescription& a, const SubtreeDescription& b )
	{
		SubtreeDescription res;
		SubtreeDescription tmp;

		for( int i = 0; i < a.SegCnt(); i++ ) {
			tmp.push_back( a[i] );
		}
		for( int i = 0; i < b.SegCnt(); i++ ) {
			tmp.push_back( b[i] );
		}
		sort( tmp.begin(), tmp.end() );
		Segment cur = tmp[0];
		for( int i = 1; i < tmp.SegCnt(); i++ ) {
			Segment u = Unite( cur, tmp[i] );
			if( u.L() == -1 ) {
				res.push_back( cur );
				cur = tmp[i];
			} else {
				swap( u, cur );
			}
		}
		if( res.empty() || cur != res.back() ) {
			res.push_back( cur );
		}
		return res;
	}
	static void BuildCondensedGraph( int cnt, const vector< int > & comp, const vector< vector< int > >& g, vector< vector< int > >& gc )
	{
		gc.assign( cnt, vector<int>() );
		for( int v = 0; v < g.size(); v++ ) {
			for( int to : g[v] ) {
				if( comp[v] == comp[to] ) {
					continue;
				}
				gc[comp[v]].push_back( comp[to] );
			}
		}
		for( int i = 0; i < gc.size(); i++ ) {
			sort( gc[i].begin(), gc[i].end() );
			auto it = std::unique( gc[i].begin(), gc[i].end() );
			gc[i].resize( std::distance( gc[i].begin(), it ) );
		}
	}
	static void GetStrongConnectivityComponents( const vector< vector< int > >& g, const vector< vector< int > >& gr, vector< int >& comp, int& cnt )
	{
		const int n = g.size();
		comp.resize( n );
		vector< int > order, tin(n), tout(n);
		TopSort( g, order, tin, tout );
		fillComponents( gr, comp, order, cnt );
	}
	static void TopSort( const vector< vector< int > >& g, vector< int >& order, vector<int>& tin, vector<int>& tout )
	{
		vector< char > used( g.size() );
		for( int i = 0; i < g.size(); i++ ) {
			if( used[i] > 0 ) {
				continue;
			}
			dfs( g, used, order, i, tin, tout);
		}
		reverse( order.begin(), order.end() );
	}

	static void Subtree( const vector< vector< int > >& g, vector< int >& order, vector<SubtreeDescription>& subTree, vector<int>& tin, vector<int>& tout)
	{
		subTree.resize( g.size() );
		for( int i = 0; i < tin.size(); i++ ) {
			subTree[i].emplace_back( tout[i], tin[i] );
		}
		vector< char > used( g.size() );
		for( int i = 0; i < g.size(); i++ ) {
			if( !used[i] )
				subtreeRecursive( g, used, i, order, subTree );
		}
	}

private:

	static void fillComponents( const vector< vector< int > >& gr, vector<int>& comp, vector< int >& order, int& cnt )
	{
		vector< char > used( gr.size() );
		cnt = 0;
		for( int i = 0; i < gr.size(); i++ ) {
			int v = order[i];
			if( used[v] > 0 ) {
				continue;
			}
			getSingleComp( gr, used, comp, v, cnt );
			cnt++;
		}
	}

	static void subtreeRecursive( const vector< vector< int > >& g, vector<char>& used,
		int v, vector<int>& order, vector<SubtreeDescription>& subTree )
	{
		// used[v] = 1;
		// for( int to : g[v] ) {
		// 	if( used[to] > 0 ) {
		// 		subTree[v] = Unite( subTree[v], subTree[to] );
		// 		continue;
		// 	}
		// 	subtreeRecursive( g, used, to, order, subTree );
		// }
		vector<int> pos(g.size());
		stack<int> st;
		st.push(v);
		while(!st.empty()){
			v = st.top();
			used[v] = 1;
			if( pos[v] < g[v].size() ){
				int to = g[v][pos[v]];
				pos[v]++;
				if( used[to] > 0 ){
					subTree[v] = Unite( subTree[v], subTree[to] );
					continue;
				}
				st.push( to );
				continue;
			}
			st.pop();
		}
	}

	static void dfs( const vector< vector< int > > & g, vector<char>& used, vector<int>& order, int v, vector<int>& tin, vector<int>& tout )
	{
		// depth++;
		// cout << depth << endl;
		// used[v] = 1;
		// tin[v] = g.size() - order.size();
		
		// for( int to : g[v] ) {
		// 	if( used[to] > 0 ) {
		// 		continue;
		// 	}
		// 	dfs( g, used, order, to, tin, tout);
		// }
		// order.push_back( v );
		// tout[v] = g.size() - order.size();
		
		vector<int> pos(g.size());
		stack<int> st;
		st.push(v);
		while(!st.empty()){
			v = st.top();
			used[v] = 1;
			if(pos[v] == 0)
				tin[v] = g.size() - order.size();
			if( pos[v] < g[v].size() ){
				int to = g[v][pos[v]];
				pos[v]++;
				if( used[to] > 0 ){
					continue;
				}
				st.push( to );
				continue;
			}
			st.pop();
			order.push_back( v );
			tout[v] = g.size() - order.size();
		}

	}

	static void getSingleComp( const vector< vector< int > > & g, vector<char>& used, vector<int>& compToVertex, int v, int compId )
	{
		// used[v] = 1;
		// for( int to : g[v] ) {
		// 	if( used[to] > 0 ) {
		// 		continue;
		// 	}
		// 	getSingleComp( g, used, compToVertex, to, compId );
		// }
		// compToVertex[v] = compId;
		int depth = 0;
		vector<int> pos(g.size());
		stack<int> st;
		st.push(v);
		while(!st.empty()){
			v = st.top();
			used[v] = 1;
			if( pos[v] < g[v].size() ){
				int to = g[v][pos[v]];
				pos[v]++;
				if( used[to] > 0 ){
					continue;
				}
				st.push( to );
				continue;
			}
			st.pop();
			compToVertex[v] = compId;
		}
	}
};


int main(int argc, char * argv[])
{
	if( argc < 3){
		cout << "use ./main input.txt output.txt";
	}
	
	const string inputFilePath( argv[1] );
	const string outputFilePath(argv[2]);
	// cin.tie( nullptr );
	// ios::sync_with_stdio( false );
// #ifdef _DEBUG
	ifstream fin( inputFilePath, std::ifstream::in) ;
	ofstream fout( outputFilePath, std::ifstream::out);
// #endif

	using namespace std::chrono;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	{


		int n = 0, m = 0;
		vector< pair<int, int> > edges;
		int x, y;
		while( fin >> x >> y ) {
			edges.emplace_back( x, y );
			n = max( n, max( x, y ) );
		}

		//cout <<  n;

		n++;
		m = edges.size();
		vector< vector< int> > g( n ), gr( n );
		vector< vector< int> > gc;
		for( int i = 0; i < m; i++ ) {
			int x, y;
			x = edges[i].first;
			y = edges[i].second;
			g[x].push_back( y );
			gr[y].push_back( x );
		}

		vector< int > compToVertex;
		vector< vector<int> > comp;
		int compCnt;
		Solution::GetStrongConnectivityComponents( g, gr, compToVertex, compCnt );
		Solution::BuildCondensedGraph( compCnt, compToVertex, g, gc );
		comp.resize( compCnt );
		for( int i = 0; i < compToVertex.size(); i++ ) {
			comp[compToVertex[i]].push_back( i );
		}
		// for( int i = 0; i < comp.size(); i++ ) {
		// 	cout << comp[i].size() << endl;
		// }

		vector< int > tsOrder, tin(gc.size()), tout(gc.size());
		Solution::TopSort( gc, tsOrder, tin, tout);
		vector< Solution::SubtreeDescription > subTree;
		//return 0;
		Solution::Subtree( gc, tsOrder, subTree, tin, tout);
				
		vector< vector<int> > blackHoles( g.size() );
		//fout << "Black holes: " << endl;
		for( int i = 0; i < subTree.size(); i++ ){
			int val = 0;
			// cout << i << ": ";
			for( int seg = 0; seg < subTree[i].SegCnt(); seg++ ){
				for( int j = subTree[i][seg].L(); j < subTree[i][seg].R(); j++ ){
					int c = tsOrder[j];
					val += comp[c].size();
				}
			}
			blackHoles[val].push_back( i );
			// cout << endl;
		}
		// for( int i =0; i < blackHoles.size(); i++ ){
		// 	if( blackHoles[i] == 0){
		// 		continue;
		// 	}
		// 	cout << i << ": " << blackHoles[i] << endl;
		// }
		// cout << endl << "gc " << endl;
		// for( int i =0; i < gc.size(); i++ ){
		// 	if( gc[i].size() == 0){
		// 		continue;
		// 	}
		// 	cout << i << ": ";
		// 	for( int x : gc[i] ){
		// 		cout << x << " ";
		// 	}
		// 	cout << endl;
		// }
		// cout << "Comp sizes " << endl;
		// for( int i = 0; i < comp.size(); i++ ){
		// 	if(comp[i].size() > 1){
		// 		cout << i << ": " << comp[i].size() << endl;
		// 	}
		// }
		// for( int i = 0 ; i < blackHoles.size(); i++ ){
		// 	cout <<  blackHoles[i].size() << endl;
		// }
		int total = 0;
		fout << inputFilePath << endl;
		fout << "Nodes " << n << "Edges " << m << endl;
		fout << "BlackHoles cnt: " << endl; 
		for( int i = 0 ; i < blackHoles.size(); i++ ){
			if( blackHoles[i].size() > 0 ){
				total += blackHoles[i].size();
				fout <<" Size: " << i << " cnt "<< blackHoles[i].size() << endl;
			}
		}
		fout << "Total: " << total << endl;
		// fout << "blackHoles: " << endl; 
		// for( int i = 0 ; i < blackHoles.size(); i++ ){
		// 	cout << i << ": " << 
		// 	for( int j =0 ; j < blackHoles[i].size(); j++ ){
		// 		int v = blackHoles[i][j];
		// 		fout << v << ": ";
		// 		for( int seg = 0; seg < subTree[v].SegCnt(); seg++ ){
		// 			for( int j = subTree[v][seg].L(); j < subTree[v][seg].R(); j++ ){
		// 				int c = tsOrder[j];
		// 				fout << c << " ";
		// 				// for( int x : comp[c] ){
		// 				// 	cout << x << " ";
		// 				// }
		// 			}
		// 		}
		// 		fout << endl;
		// 	}
		// }

		//cout << "TsOrder" << endl;
		//for( int x : tsOrder ) {
		//	cout << x << " ";
		//}
		//cout << endl;
		//for( int i = 0; i < subTree.size(); i++ ) {
		//	cout << i << ": ";
		//	auto x = subTree[i];
		//	for( auto xx : x ) {
		//		cout << xx.L() << " " << xx.R() << "; ";
		//	}
		//	cout << endl;
		//}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>( t2 - t1 );
	fout << "It took me: "<< time_span.count() <<  " seconds" << endl;
	cout << "It took me: "<< time_span.count() <<  " seconds" << endl;
	return 0;
}