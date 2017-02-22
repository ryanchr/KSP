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

#include <fstream>
#include <cstdlib>

using namespace std;





int main(int argc, char * argv[])
{
    int row_num = 60;
    int column_num = 5; // tmp fix it as 5 in the grid network since the max_phy_hops is 7.
    int total_node_num = row_num * column_num;


    ofstream myfile ("randomGridData/vertex.csv");
    if (myfile.is_open())
    {
        myfile << "id|nodeName\n";
        for (int idx = 0; idx < total_node_num; ++idx)
        {
            myfile << idx << "|noName\n";
        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open file";
    } 
        

    ofstream myfile_edge ("randomGridData/edge.csv");
    if (myfile_edge.is_open())
    {
        myfile_edge << "Node.id|Node.id|Cost|Distance\n";

        for (int idx = 0; idx < total_node_num; ++idx)
        {
            float costValue;
            float distanceValue;
            int dest_id;
            float costRange = 1; // rand float between 0~10
            float distanceRange = 0.1; // rand float between 0~10

            if (idx%60 == 0)
            {
                int dest_id_array[5];
                dest_id_array[0] = idx - row_num;
                dest_id_array[1] = idx - row_num + 1;
                dest_id_array[2] = idx + 1;
                dest_id_array[3] = idx + row_num;
                dest_id_array[4] = idx + row_num + 1;
                for (int idx2 = 0; idx2 < 5; ++idx2)
                {
                    dest_id = dest_id_array[idx2];
                    if (dest_id>=0 && dest_id<total_node_num)
                    {
                        myfile_edge << idx << "|" << dest_id<< "|";
                        costValue     = rand()/(float)RAND_MAX * costRange;   // rand float between 0~10
                        distanceValue = rand()/(float)RAND_MAX * distanceRange;   // rand float between 0~3 
                        myfile_edge << costValue << "|" << distanceValue << "\n";
                    }                
                }
            }
            else if (idx%60 == 59)
            {
                int dest_id_array[5];
                dest_id_array[0] = idx - row_num - 1;
                dest_id_array[1] = idx - row_num;
                dest_id_array[2] = idx - 1;
                dest_id_array[3] = idx + row_num - 1;
                dest_id_array[4] = idx + row_num;
                for (int idx2 = 0; idx2 < 5; ++idx2)
                {
                    dest_id = dest_id_array[idx2];
                    if (dest_id>=0 && dest_id<total_node_num)
                    {
                        myfile_edge << idx << "|" << dest_id<< "|";
                        costValue     = rand()/(float)RAND_MAX * costRange;   // rand float between 0~10
                        distanceValue = rand()/(float)RAND_MAX * distanceRange;   // rand float between 0~3 
                        myfile_edge << costValue << "|" << distanceValue << "\n";
                    }                
                }
            }
            else
            {
                int dest_id_array[8];
                dest_id_array[0] = idx - row_num - 1;
                dest_id_array[1] = idx - row_num;
                dest_id_array[2] = idx - row_num + 1;
                dest_id_array[3] = idx - 1;
                dest_id_array[4] = idx + 1;
                dest_id_array[5] = idx + row_num - 1;
                dest_id_array[6] = idx + row_num;
                dest_id_array[7] = idx + row_num + 1;
                for (int idx2 = 0; idx2 < 8; ++idx2)
                {
                    dest_id = dest_id_array[idx2];
                    if (dest_id>=0 && dest_id<total_node_num)
                    {
                        myfile_edge << idx << "|" << dest_id<< "|";
                        costValue     = rand()/(float)RAND_MAX * costRange;   // rand float between 0~10
                        distanceValue = rand()/(float)RAND_MAX * distanceRange;   // rand float between 0~3 
                        myfile_edge << costValue << "|" << distanceValue << "\n";
                    }                
                }
            }

        }

        myfile_edge.close();
    }
    else
    {
        cout << "Unable to open file";
    } 

    return 0;
}  // end main

