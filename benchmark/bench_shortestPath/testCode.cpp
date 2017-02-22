//====== Graph Benchmark Suites ======//
//========== Shortest Path ===========//
// 
// Single-source shortest path
// 
// Usage: ./sssp    --dataset <dataset path> 
//                  --root <root vertex id> 
//                  --target <target vertex id>

#include <cassert>
#include <iostream>

#include <vector>
#include <list>
#include <string>
#include <queue> // for queue and priority_queue 
#include <stack>
#include <map> //for map and multimap
#include <unordered_map> //for unordered_map and unordered_multimap
#include <algorithm>
using namespace std;




void displayVectorCont(const vector< pair<int,float> > & vecInput)
{
    for (auto idx = vecInput.begin(); idx != vecInput.end(); ++idx)
        cout << idx->first << "->" << idx->second << endl; 
}

int main(int argc, char * argv[])
{

    vector< pair<int,float> > vector_of_pair; 
    vector_of_pair.push_back(pair<int,float>(4,6.3));
    vector_of_pair.push_back(pair<int,float>(2,9.5));
    vector_of_pair.push_back(pair<int,float>(1,19.5));
    vector_of_pair.push_back(pair<int,float>(3,12.5));

    displayVectorCont(vector_of_pair);
    cout<<"after sort\n";
/*    sort(vector_of_pair.begin(),vector_of_pair.end(),
        [](pair<int,float> vecInput_1, pair<int,float> vecInput_2) -> bool
        {
            return (vecInput_1.second < vecInput_2.second);
        }
        );*/
    sort(vector_of_pair.begin(),vector_of_pair.end(),
        [](pair<int,float> vecInput_1, pair<int,float> vecInput_2) {return (vecInput_1.second < vecInput_2.second);} );
        
            
        
            
    displayVectorCont(vector_of_pair);
    cout<<"output based on index"<<endl;
    for (int idx=0;idx<vector_of_pair.size();++idx)
        cout << vector_of_pair[idx].first << "->" << vector_of_pair[idx].second << endl; 
    cout<<"==================================================================\n";

    #define MY_INFINITY 0xfff0
    uint64_t dd = MY_INFINITY;
    if (dd==MY_INFINITY)
        cout<<"asdfasfasdfas"<<endl;
    string sss="25";
    int s2ivalue = atoi(sss.c_str());
    cout<<s2ivalue<<endl;
    cout<<"==================================================================\n";
    assert(true);

    return 0;
}  // end main

