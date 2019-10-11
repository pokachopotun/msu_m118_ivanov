#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
using namespace std;

const float ThermoSpeed = 0.25; 
const float Cooling = 0.95;
const float AntsHeating = 0.1;

class Node {
public:
	int Id;
	float Heat;
	float NewHeat;
	int AntsCnt;
	int RecvAntsCnt;
	vector<Node*> Succ;
	vector<Node*> Pred;
	Node( int idx = -1 ) : Id(idx), Heat(0), NewHeat(0), AntsCnt(0), RecvAntsCnt(0), Succ(), Pred()
	{
	}
	
	void CalcNewHeat_old() {
		float sumHeat = Heat;
		int cnt = 1;
		
		for (Node* node : Succ) {
			const float heat = node->Heat;
			if (heat < Heat){
				sumHeat += heat;
				cnt += 1;
			}
		}
	
		for (Node* node : Pred) {
			const float heat = node->Heat;
			if (heat < Heat){
				sumHeat += heat;
				cnt += 1;
		}
		}

		const float targetHeat = sumHeat / cnt;
		NewHeat = Heat + ( targetHeat - Heat ) * ThermoSpeed;
	}

	void CalcNewHeat() {
		float sumHeat = Heat;
		int cnt = 1;
		
		for (Node* node : Succ) {
			const float heat = node->Heat;
			if (heat < Heat){
				sumHeat += heat;
				cnt += 1;
			}
		}
	
		for (Node* node : Pred) {
			const float heat = node->Heat;
			if (heat < Heat){
				sumHeat += heat;
				cnt += 1;
			}
		}

		const float targetHeat = sumHeat / cnt;
		NewHeat = Heat + ( targetHeat - Heat ) * ThermoSpeed;
	}
	
	void Cool(){
		Heat *= Cooling;	
	}
	
	void HeatAnts(){
		Heat += AntsCnt * AntsHeating;
	}	

	void UpdateHeat() {
		Heat = NewHeat;
	}
	
	void Spawn(int cnt) {
		AntsCnt += cnt;
	}
	
	void SendAnts(){
		float maxHeat = 0;
		int cnt = 0;
		if (Succ.size() == 0)
			return;
		for( Node* node : Succ ){	
			maxHeat = max( maxHeat, node->Heat );
		}

		vector<float> p;
		float psum = 0;
		for( Node* node : Succ) {
			p.push_back(maxHeat - node->Heat);
			psum += p.back();
		}
		if (psum == 0.0) {
			for (float& pp : p) {
				pp = 1.0 / p.size();
			}
		} else {
			for (float& pp : p) {
				pp /= psum;
			}
		}

		vector<float> intervals(p.size());
		for( int i = 0; i < intervals.size(); i++ ) {
			intervals[i] = i;
		}
		std::default_random_engine generator;
		std::piecewise_linear_distribution<float> distribution(intervals.begin(), intervals.end(), p.begin());
		
		int antsSent = 0;	
		for (int i = 0 ; i < AntsCnt; i++ ) {
			int pos = distribution(generator);
			if (Succ[pos]->Heat <= Heat) {
			//	#pragma omp critical 
				{
					Succ[pos]->RecvAntsCnt++;
				}
				antsSent++;
			}	
		}
		AntsCnt -= antsSent;
	}

	void RecvAnts() {
		AntsCnt += RecvAntsCnt;
		RecvAntsCnt = 0;
	}

	void PrintNode() {
		cout << "Node Id " << Id << endl;
		cout << "Succ: ";
		for(Node* node : Succ) {
			cout << node->Id << " ";
		}
		cout << endl << "Pred: ";
		for(Node* node : Pred) {
			cout << node->Id << " ";
		}
		cout << endl;
	}
};


class Solution {
public:
	enum ESource {
		S_None = -1,
		S_File_Text = 0,
		S_File_Binary = 1
	};

	int Iter;
	vector<Node*> Nodes;
	
	Solution(int nThreads = 1) : Iter(0), Nodes()
	{
	}

	~Solution() {
		for (Node* node : Nodes) {
			if (node != nullptr) {
				delete node;
			}
		}
	}

	void ReadGraphFromTextFile(const string& graphFilePath) {
		ifstream ifs(graphFilePath);
		int nNodes, nLines;
		ifs >> nNodes >> nLines;
		for (Node * node : Nodes){
			if (node != nullptr)
				delete node;
		}

		Nodes.assign(nNodes, nullptr);
		
		for (int i =0 ; i < Nodes.size(); i++ ) {
			Nodes[i] = new Node(i);
		}
		
		for (int line = 0; line < nLines; line++ ) {
			int f, n, a;
			ifs >> f >> a >> n;
			Nodes[f]->AntsCnt = a;
			for (int i = 0; i < n; i++) {
				int to;
				ifs >> to;
				Nodes[f]->Succ.push_back(Nodes[to]);
				Nodes[to]->Pred.push_back(Nodes[f]);
			}
		}
	}

	void Diffuse() {
		
		for(int i = 0; i < Nodes.size(); i++) {
			Nodes[i]->CalcNewHeat();
		}
		for(int i = 0; i < Nodes.size(); i++) {
			Nodes[i]->UpdateHeat();
		}
	}

	void Cool(){
		for(int i = 0; i < Nodes.size(); i++) {
			Nodes[i]->Cool();
		}
	}

	void MoveAnts() {
		for(int i = 0; i < Nodes.size(); i++) {
			Nodes[i]->SendAnts();
		}
		for(int i = 0; i < Nodes.size(); i++) {
			Nodes[i]->RecvAnts();
		}
	}
	
	void Heat(){
		for(int i = 0; i < Nodes.size(); i++) {
			Nodes[i]->HeatAnts();
		}
	}
	
	void RunStep() {
		MoveAnts();
		Heat();
		Diffuse();
		Cool();
//		PrintState();
	}

	void Run(int numIter){
//		CheatInit();
	//	PrintState();
		while (numIter--) {
			RunStep();
	//		PrintState();
//			Show();
//			SaveFigure();
		}
	}

	void Show() {
		cerr << " Solution::Show() Not implemented" << endl;
	}

	void SaveFigure() {
		cerr << " Solution::SaveFigure() not implementeed" << endl;
	}

	void PrintGraph() {
		for (Node* node : Nodes) {
			node->PrintNode();
		}
	}

	void PrintState() {
		cout << "Iter " << Iter;
		for (Node* node : Nodes) {
			cout << " Node " << node->Id;
			cout << " Heat " << node->Heat;
			cout << " AntsCnt " << node->AntsCnt;
			cout << endl;
		} 
	}
	
	void CheatInit() {
		for ( Node* node : Nodes ) {
			node->Heat = 10.0 / ( node->Id + 1) ;
		}
	}
};

int main(int argc, char* argv[]){
	if (argc < 4){
		cout << "Use ./main graphFilePath nIter nThreads" << endl;
		return 0;
	}
	const string graphFilePath(argv[1]);
	const int nIter(atoi(argv[2]));
	const int nThreads(atoi(argv[3]));
	omp_set_num_threads(nThreads);
	Solution sln;
	sln.ReadGraphFromTextFile(graphFilePath);
	//sln.//PrintGraph();
	double time = omp_get_wtime();
	sln.Run(nIter);
	time = omp_get_wtime() - time;
	int tid = omp_get_thread_num();
	if ( tid == 0 ) {
		cout << "graph " <<  graphFilePath << " nIter " << nIter << " nThreads " << nThreads << " time " << time << endl;
	}
	return 0;
}
