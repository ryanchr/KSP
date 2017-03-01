//====== Graph Benchmark Suites ======//
//========== Shortest Path ===========//
// 
// Single-source shortest path
// 
// Usage: ./sssp    --dataset <dataset path> 
//                  --root <root vertex id> 
//                  --target <target vertex id>

#include "common.h"
#include "def.h"
#include "openG.h"
#include <queue>
#include "omp.h"
#include <fstream>

#define OUTPUT_PATHS

//#define TOPKSP_DEBUG
#define TOPKSP_PRINTOUT


#ifdef HMC
#include "HMC.h"
#endif

#ifdef SIM
#include "SIM.h"
#endif

#define MY_INFINITY 0xfff0  // for unsigned i.e., 65520

#define ENABLE_OUTPUT  //new
#include <cfloat>      //new  for FLT_MAX
#include <cassert>
#include <stack>
#include <unordered_set>
#include <cstdlib>

using namespace std;
size_t beginiter = 0;
size_t enditer = 0;

typedef pair<size_t,float> pair_IntFlt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<float,size_t> pair_FltInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<size_t,size_t> pair_IntInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
/*class vertex_property
{
public:
    vertex_property():min_cost(MY_INFINITY),predecessor(MY_INFINITY),update(MY_INFINITY){}

    uint16_t min_cost;
    uint64_t predecessor;
    uint16_t update;
};*/
/*class edge_property
{
public:
    edge_property():cost(1){}

    uint16_t cost;

};*/
class vertex_property // new
{
public:  // sum_distance(0),sum_hops(0){}
    vertex_property():min_cost(FLT_MAX),successor(MY_INFINITY),sum_distance(FLT_MAX),sum_hops(MY_INFINITY){}

    float min_cost;    // for shortest path
    // predecessor; // for shortest path successor,
    uint64_t successor; // for shortest path 

    float sum_distance;    // new
    uint64_t sum_hops;    // new 



    //uint64_t update;

    vector<pair_IntFlt> sorted_edges_of_vertex;  // pair: id (outedge), reduced_cost
};
class edge_property // new
{
public:
    edge_property():cost(FLT_MAX),phy_dist(FLT_MAX),reduced_cost(FLT_MAX){} // new:  Note: 

    float cost; // note: here cost means cost
    float phy_dist; //new

    float reduced_cost; //  
};
typedef openG::extGraph<vertex_property, edge_property> graph_t;
typedef graph_t::vertex_iterator    vertex_iterator;
typedef graph_t::edge_iterator      edge_iterator;



void reset_graph(graph_t & g)
{
    vertex_iterator vit;
    for (vit=g.vertices_begin(); vit!=g.vertices_end(); vit++)
    {
        vit->property().min_cost     = FLT_MAX;
        vit->property().successor    = MY_INFINITY;
        vit->property().sum_distance = FLT_MAX;
        vit->property().sum_hops     = MY_INFINITY;
        vit->property().sorted_edges_of_vertex.clear();

        for (edge_iterator eit = vit->in_edges_begin(); eit != vit->in_edges_end(); eit++)  // new for in edge
        {
            eit->property().reduced_cost = FLT_MAX;
        }
    }
}







class vertex_property_tau 
{
public:   
    vertex_property_tau():at_KSPaths(false),predecessor(MY_INFINITY),internel_id_of_g(MY_INFINITY),min_cost(0),sum_distance(0),sum_hops(0){}

    bool at_KSPaths; //vector<size_t> KSPaths_record;
    uint64_t predecessor; // for shortest path, store the id of the graph tau
    uint64_t internel_id_of_g; // internel id of orginal graph g

    float min_cost;    // for shortest path
    float sum_distance;    // for shortest path
    size_t sum_hops;    // for shortest path
};
class edge_property_tau // new
{
public:
    edge_property_tau(){} //edge_property_tau():cost(FLT_MAX),phy_dist(FLT_MAX),reduced_cost(FLT_MAX){} 
    //float cost; 
    //float phy_dist;
};
typedef openG::extGraph<vertex_property_tau, edge_property_tau> graph_tau;
typedef graph_tau::vertex_iterator    vertex_iterator_tau;
typedef graph_tau::edge_iterator      edge_iterator_tau;







/*class comp
{
public:
    bool operator()(pair_IntFlt a, pair_IntFlt b)
    {
        return a.second > b.second;
    }
};*/
class min_comp_IntFlt
{
public:
    bool operator()(pair_IntFlt a, pair_IntFlt b)
    {
        return a.second > b.second;
    }
};
class min_comp_FltInt
{
public:
    bool operator()(pair_FltInt a, pair_FltInt b)
    {
        return a.first > b.first;
    }
};


#ifdef TOPKSP_DEBUG
    template <typename T>
    void displayVectorOrListContents(const T& input)
    {
        for (auto idx = input.cbegin(); idx != input.cend(); ++ idx)
        {
            cout << *idx << ",";
        }
        cout << endl;
    }
    template <typename T>
    void displayVectorOrListOfPairContents(const T& input)  
    {
        for (auto idx = input.cbegin(); idx != input.cend(); ++ idx)
        {
            cout << "(" << idx->first << "," << idx->second << ")" << ",";
        }
        cout << endl;
    }


    void output_shorest_path(graph_t& g, size_t src, size_t dest) //for test, src and dest is exID, dest is the input dest of top_ksp func.
    {
        cout<<"the shortest path is: ";
        uint64_t internel_dest = g.external_to_internel_id(to_string(dest)); 
        uint64_t curr_id_g = g.external_to_internel_id(to_string(src)); 

        uint64_t curr_exID = src;
        cout<<curr_exID;
        do
        {
            vertex_iterator vit = g.find_vertex(curr_id_g);   
            curr_id_g  = vit->property().successor; 

            if (curr_exID==289)
            {
                cout<<"curr_id_g="<<curr_id_g<<endl;
            }

            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            cout<<"-->"<<curr_exID;
            //cout<<"-->"<<curr_exID<<endl;
        } while (curr_id_g != internel_dest);
        cout<<endl;
    }

/*    void output_shorest_path(graph_t& g, size_t src, size_t dest) //for test, src and dest is exID, dest is the input dest of top_ksp func.
    {
        cout<<"the shortest path is: ";

        uint64_t internel_dest = g.external_to_internel_id(to_string(dest)); 
        uint64_t curr_id_g = g.external_to_internel_id(to_string(src)); 

        cout << "curr_id_g = " << curr_id_g <<endl;


        uint64_t curr_exID = src;
        cout<<curr_exID;
        do
        {
            vertex_iterator vit = g.find_vertex(curr_id_g); 
            curr_id_g  = vit->property().successor; 

            cout<<"after vit->property().successor at curr_id_g = "<<curr_id_g<<endl;
            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            cout<<"after atoi at new curr_exID = "<<curr_exID<<endl;

            //cout<<"-->"<<curr_exID;
            
        } while (curr_id_g != internel_dest);
        cout<<endl;
    }*/
#endif


vertex_iterator_tau add_partialSP_totau(graph_t& g, size_t src_id_g, size_t dest_id_g, graph_tau& tau, size_t start_id_tau, double max_phy_dist, size_t max_phy_hops)  
{
    vertex_iterator_tau src_vit_tau;
    vertex_iterator_tau dest_vit_tau;
    edge_iterator_tau eit_tau;
    vertex_iterator src_vit_g;
    vertex_iterator dest_vit_g;
    edge_iterator eit_g;

    src_vit_tau = tau.find_vertex(start_id_tau);  
    assert(src_vit_tau->property().internel_id_of_g == src_id_g); // this is a require for the following code to work

    uint64_t tau_tmp_id = start_id_tau;
    src_vit_g = g.find_vertex(src_id_g);  
    while (src_vit_g->id() != dest_id_g) // note: dest_id_g should be the internel_dest
    {
        dest_vit_g = g.find_vertex( src_vit_g->property().successor ); // new, ori is  predecessor


        //cout << "src_vit_g->id()= "  << src_vit_g->id()  << endl;
        //cout << "dest_vit_g->id()= " << dest_vit_g->id() << endl;
        bool find_result = g.find_out_edge_2id(src_vit_g->id(), dest_vit_g->id(), eit_g); // for eit_g->property().cost and eit_g->property().phy_dist
        assert(find_result);



        // add point and edge at tau
        dest_vit_tau = tau.add_vertex();
        dest_vit_tau->property().at_KSPaths = false;// 
        dest_vit_tau->property().predecessor = tau_tmp_id; // 
        dest_vit_tau->property().internel_id_of_g = dest_vit_g->id();   
        dest_vit_tau->property().min_cost     = src_vit_tau->property().min_cost + eit_g->property().cost; //   
        dest_vit_tau->property().sum_distance = src_vit_tau->property().sum_distance + eit_g->property().phy_dist; //   
        dest_vit_tau->property().sum_hops     = src_vit_tau->property().sum_hops + 1;

        //cout<<"the new added node (idx_g,idx_tau) in Tau is ("<<dest_vit_tau->property().internel_id_of_g<<","<<dest_vit_tau->id()<<")"<<endl;
        //cout<<"The new added node (min_cost,sum_distance,sum_hops) in Tau is ("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";
        //cout<<endl;

        tau.add_edge(src_vit_tau->id(),dest_vit_tau->id(),eit_tau); // note: put all info at vertex, see dest_vit_tau
        
        //new   
        ////if (dest_vit_tau->property().sum_distance >= max_phy_dist || dest_vit_tau->property().sum_hops >= max_phy_hops)
        //// this following is only for single-layer real data
        if (dest_vit_tau->property().sum_hops >= max_phy_hops)
            break;

        // for next iteration use
        tau_tmp_id = dest_vit_tau->id();  
        src_vit_g  = dest_vit_g; 
        src_vit_tau  = dest_vit_tau;  
    }
    return dest_vit_tau;
}
 


