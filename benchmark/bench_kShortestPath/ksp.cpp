//====== Graph Benchmark Suites ======//
//========== Shortest Path ===========//
// 
// Single-source shortest path
// 
// Usage: ./sssp    --dataset <dataset path> 
//                  --root <root vertex id> 
//                  --target <target vertex id>

#include "omp.h"
#include <queue>
#include <cfloat>     
#include <cassert>
#include <stack>
#include <cstdlib>
using namespace std;

vector<vector<int> > ksp(vector<pair<int, int> >edges, pair<int,int> num_hops_range  \
					   , int src, vector<int> dest, int num_virtexs)
{
	//Initialization
	unordered_map<int, int> dis, pre;
	vector<vector int> res;
	dis[src] = 0;

	//
	

}