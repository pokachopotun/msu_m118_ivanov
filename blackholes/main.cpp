#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <queue>

using namespace std;

using graph_t = vector< vector< int > * >;
using blackhole_t = vector< int >;

graph_t parseGraphFile(const string& filePath){
	graph_t g, g_rev;
	ifstream gin(filePath);
	int n, m;
	gin >> n >> m;
	for(int i = 0; i < n; i++){
		g.push_back(new vector<int>());
		g_rev.push_back(new vector<int>());
	}
	for(int i= 0; i < m; i++){
		int x,y;
		cin >> x >> y;
		g[x]->push_back(y);
		g_rev[y]->push_back(x);
	}
	return g;
}

vector<int> getSubtree(const graph_t& g, int s){
	vector<int> ans;
	queue<int> q;
	q.push(s);
	vector< char > used(g.size(), 0);
	used[s] = 1;
	while(!q.empty()){
		int v = q.front();
		q.pop();
		ans.push_back(v);
		for(int to : *g[v]){
			if(used[to] || g[to] == nullptr ){
				continue;
			}
			used[to] = true;
			q.push(to);
		}
	}
	return ans;
}

void dye(const graph_t& g, int s, vector<int>& used, int color){
	queue<int> q;
	q.push(s);
	used[s] = color;
	while(!q.empty()){
		int v = q.front();
		q.pop();
		for(int to : *g[v]){
			if( g[to] == nullptr ){
				continue;
			}
			if(used[to] > 0){
				used[to];
			}
			used[to] = color;
			q.push(to);
		}
	}
}

graph_t* getSubgraph(const graph_t& g, int n){
	graph_t* s = new graph_t;
	s->assign(g.size(), nullptr);
	for(int i = 0; i < g.size(); i++){
		if( g[i]->size() < n ){
			(*s)[i] = g[i];
		}
	}
	return s;
}

graph_t* getReverseSubGraph(const graph_t& gRev, const graph_t& sub){
	graph_t* res = new graph_t( gRev );
	for(int v = 0; v < gRev.size(); v++){
		if( sub[v] == nullptr ){
			(*res)[v] = nullptr;
		}
	} 
	return res;
}

vector< blackhole_t > findBlackHoles(const graph_t& g, int n){
	n++;
	vector< blackhole_t > ans;
	graph_t* potential = nullptr, *potentialPrev = nullptr, *potentialRev = nullptr;
	for( int bHoleSize = 1; bHoleSize < n; bHoleSize++){
		

		// build potential list
		delete potentialPrev;
		potentialPrev = potential;
		potential = getSubgraph(g, bHoleSize);
		potentialRev = getReverseSubGraph(g, *potential);

		//select candidate list
		for(int v = 0; v < potential->size(); v++){
			if((*potential)[v] == nullptr){
				continue;
			}
			for(int to : *( (*potential)[v] ) ){
				if((*potentialPrev)[v] != nullptr && (*potential)[to] == nullptr){
					vector< int >  subTree = getSubtree(*potentialRev, v);
					for(int rm : subTree){
						(*potential)[rm] = (*potentialRev)[rm] = nullptr;
					}
				}
			}
		}

		//select candidate list
		for(int v = 0; v < potential->size(); v++){
			if((*potential)[v] == nullptr){
				continue;
			}
			vector< int >  subTree = getSubtree(*potentialRev, v);
			if( subTree.size() > bHoleSize ){
				for(int rm : subTree){
					(*potential)[rm] = (*potentialRev)[rm] = nullptr;
				}
			}
			if( subTree.size() == bHoleSize ){
				ans.push_back(subTree);
				for(int rm : subTree){
					(*potential)[rm] = (*potentialRev)[rm] = nullptr;
				}
			}
		}
	}
	return ans;
}	
int main(int argc, char* argv[]){
	if(argc < 3){
		cout << "use ./blackholes input.txt output.txt" << endl;
		return 0;
	}
	const string& inputFilePath(argv[1]);
	const string& outputFilePath(argv[2]);


	return 0;
}