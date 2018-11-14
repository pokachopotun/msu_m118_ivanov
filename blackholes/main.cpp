#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <queue>

using namespace std;

using graph_t = vector< vector< int > * >;
using blackhole_t = vector< int >;

pair< graph_t, graph_t > parseGraphFile(const string& filePath){
	graph_t g, gRev;
	ifstream gin(filePath);
	int n, m;
	gin >> n >> m;
	for(int i = 0; i < n; i++){
		g.push_back(new vector<int>());
		gRev.push_back(new vector<int>());
	}
	for(int i= 0; i < m; i++){
		int x,y;
		gin >> x >> y;
		g[x]->push_back(y);
		gRev[y]->push_back(x);
	}
	return make_pair(g, gRev);
}

vector<int> getSubtree(const graph_t& g, int s){
	vector<int> ans;
	queue<int> q;
	q.push(s);
	vector< bool > used(g.size());
	used[s] = true;
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

void weakConnectivity(int s, const graph_t& g, const graph_t& gRev, vector<int>& used){
	
	queue<int> q;
	q.push(s);
	used[s] = s;
	while(!q.empty()){
		int v = q.front();
		q.pop();
		for(int to : *g[v]){
			if( g[to] == nullptr ){
				continue;
			}
			if(used[to] >= 0){
				continue;
			}
			used[to] = s;
			q.push(to);
		}
		for(int to : *gRev[v]){
			if( gRev[to] == nullptr ){
				continue;
			}
			if(used[to] >= 0 ){
				continue;
			}
			used[to] = s;
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

vector< blackhole_t > findBlackHoles(const graph_t& g, const graph_t& gRev, int n){
	n++;
	vector< blackhole_t > ans;
	graph_t* potential = nullptr, *potentialPrev = nullptr, *potentialRev = nullptr;
	for( int bHoleSize = 1; bHoleSize < n; bHoleSize++){
		cout << bHoleSize << endl;

		// build potential list
		potentialPrev = potential;
		potential = getSubgraph(g, bHoleSize);
		potentialRev = getReverseSubGraph(gRev, *potential);

		//select candidate list
		for(int v = 0; v < potential->size(); v++){
			if((*potential)[v] == nullptr){
				continue;
			}
			for(int to : *( (*potential)[v] ) ){
				if( (potentialPrev == nullptr || 
					( potentialPrev != nullptr && (*potentialPrev)[v] == nullptr )) && (*potential)[to] == nullptr){
					vector< int >  subTree = getSubtree(*potentialRev, v);
					for(int rm : subTree){
						(*potential)[rm] = (*potentialRev)[rm] = nullptr;
					}
				}
			}
		}

		//select final candidate list
		for(int v = 0; v < potential->size(); v++){
			if((*potential)[v] == nullptr){
				continue;
			}
			vector< int >  subTree = getSubtree(*potential, v);
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
		vector<int> used(g.size(), -1);
		for(int v = 0; v < potential->size(); v++){
			if((*potential)[v] == nullptr){
				continue;
			}
			if(used[v] == -1){
				weakConnectivity(v, *potential, *potentialRev, used);
			}
		}
		cout << bHoleSize << " Used:" << endl;
		for(int i=0; i < used.size(); i++){
			cout << i << " " << used[i] << endl;
		}
		cout << "<<<<<<<<<<<" << endl;
		for(int i = 0; i < potential->size(); i++){
			vector< int > comp;
			for(int v = 0; v < potential->size(); v++){
				if(used[v] == i){
					comp.push_back(v);
				}
			}
			// if(comp.size() == bHoleSize){
				ans.push_back(comp);
			// }
		}
	}
	return ans;
}

void printOutput(const vector< blackhole_t >& ans, const string& outputFilePath){
	ofstream out(outputFilePath);
	for(auto vv : ans){
		out << vv.size() << " ";
		for(auto v : vv){
			out << v << " ";
		}
		out << endl;
	}
}
int main(int argc, char* argv[]){
	if(argc < 3){
		cout << "use ./blackholes input.txt output.txt" << endl;
		return 0;
	}
	const string& inputFilePath(argv[1]);
	const string& outputFilePath(argv[2]);

	auto g = parseGraphFile(inputFilePath);
	cout << "Parsing done" << endl;
	auto ans = findBlackHoles(g.first, g.second, g.first.size());

	printOutput(ans, outputFilePath);
	return 0;
}