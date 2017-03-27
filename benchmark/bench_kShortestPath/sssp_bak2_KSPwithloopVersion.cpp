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

//#define TOPKSP_DEBUG



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

using namespace std;
size_t beginiter = 0;
size_t enditer = 0;

typedef pair<size_t,float> pair_IntFlt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<float,size_t> pair_FltInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //

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
public:
    vertex_property():min_cost(FLT_MAX),successor(MY_INFINITY){}

    float min_cost;    // for shortest path
    // predecessor; // for shortest path successor
    uint64_t successor; // for shortest path 

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



class vertex_property_tau 
{
public:   
    vertex_property_tau():at_KSPaths(false),predecessor(MY_INFINITY),internel_id_of_g(MY_INFINITY),min_cost(0),sum_distance(0),sum_hops(0){}

    bool at_KSPaths; //vector<size_t> KSPaths_record;
    uint64_t predecessor; // for shortest path, store the id of the graph tau
    uint64_t internel_id_of_g; // internel id of orginal graph g

    float min_cost;    // for shortest path
    float sum_distance;    // for shortest path
    uint64_t sum_hops;    // for shortest path
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
            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            cout<<"-->"<<curr_exID;
        } while (curr_id_g != internel_dest);
        cout<<endl;
    }
#endif


vertex_iterator_tau add_partialSP_totau(graph_t& g, size_t src_id_g, size_t dest_id_g, graph_tau& tau, size_t start_id_tau)  
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

        // for next iteration use
        tau_tmp_id = dest_vit_tau->id();  
        src_vit_g  = dest_vit_g; 
        src_vit_tau  = dest_vit_tau;  
    }
    return dest_vit_tau;
}

void top_ksp(graph_t& g, size_t src, size_t dest, size_t Kvalue, gBenchPerf_event & perf, int perf_group) // src and dest are exID
{
    perf.open(perf_group);
    perf.start(perf_group);
#ifdef SIM
    SIM_BEGIN(true);
#endif
    uint64_t internel_src  = g.external_to_internel_id(to_string(src)); 
    uint64_t internel_dest = g.external_to_internel_id(to_string(dest)); 

    //// Now all processing is based on internel id first.
    //// (1) sssp sub-procedure  -->for vertex: update  v_vit->property().min_cost and v_vit->property().predecessor 
    priority_queue<pair_IntFlt, vector<pair_IntFlt>, min_comp_IntFlt> PQ;

    vertex_iterator dest_vit = g.find_vertex(internel_dest);  // note: here the source of sssp is internel_dest, 
    dest_vit->property().min_cost = 0;
    PQ.push(pair_IntFlt(internel_dest,0));
    // vit->property().predecessor is used to construct sssp, where the ancestor of all nodes is internel_dest, 
    // by using vit->property().predecessor recursively, we can find the shortest path of any node to internel_dest
/*    while (!PQ.empty()) 
    {
        size_t u = PQ.top().first; //id
        PQ.pop();

        vertex_iterator u_vit = g.find_vertex(u);
        for (edge_iterator eit = u_vit->edges_begin(); eit != u_vit->edges_end(); eit++)  // for out edge
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);
            // for every  vertex, try relaxing the path
            float alt = u_vit->property().min_cost + eit->property().cost; // 
            if (alt < v_vit->property().min_cost) 
            {
                v_vit->property().min_cost = alt; 
                v_vit->property().predecessor = u;
                PQ.push(pair_IntFlt(v,alt));
            }
        }
    }*/
/*    //debug: passed : 3-->2-->4-->5-->6-->10-->11-->12
    {
        cout<<"the shortest path is (reversed): "<<endl;
        uint64_t curr_id_g = internel_src;
        uint64_t curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
        cout<<curr_exID;
        do
        {
            vertex_iterator vit = g.find_vertex(curr_id_g);   
            curr_id_g  = vit->property().predecessor; 
            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            cout<<"-->"<<curr_exID;
        } while (curr_id_g != internel_dest);
    }
    cout<<endl;*/


    while (!PQ.empty()) 
    {
        size_t u = PQ.top().first; //id
        PQ.pop();

        vertex_iterator u_vit = g.find_vertex(u);
        for (edge_iterator eit = u_vit->in_edges_begin(); eit != u_vit->in_edges_end(); eit++)  // new for in edge
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);
            // for every  vertex, try relaxing the path
            float alt = u_vit->property().min_cost + eit->property().cost; // 
            if (alt < v_vit->property().min_cost) 
            {
                v_vit->property().min_cost = alt; 
                v_vit->property().successor = u;  // new, ori is predecessor
                PQ.push(pair_IntFlt(v,alt));
            }
        }
    }

#ifdef TOPKSP_DEBUG
    // output externel ID
    /*for (int src_tmp=0; src_tmp < g.num_vertices()-1; src_tmp++)//for MPS alg simple ex. note: cannot output dest to dest
    {
        output_shorest_path(g, src_tmp, dest);
    }//may terminate if there is no path between to nodes
    */
    // output internel ID
/*    for (int src_tmp=0; src_tmp < g.num_vertices()-1; src_tmp++)  
    {
        cout<<"the shortest path is (internel): ";
        uint64_t curr_id_g = src_tmp;
        uint64_t curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
        cout<<curr_exID;
        do
        {
            vertex_iterator vit = g.find_vertex(curr_id_g);   
            curr_id_g  = vit->property().successor; 
            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            cout<<"-->"<<curr_exID;
        } while (curr_id_g != internel_dest);
        cout<<endl;
    }*/
 #endif



    
    ////(2) reduced_cost computing procedure  -->for edge: update eit->property().reduced_cost
    ////(3) rearrange the arcs
    for (vertex_iterator u_vit=g.vertices_begin(); u_vit!=g.vertices_end(); u_vit++) // for each vertex u
    {
        for (edge_iterator eit = u_vit->edges_begin(); eit != u_vit->edges_end(); eit++) // for each outedge u->v from vertex u
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);
            float reduced_cost = v_vit->property().min_cost - u_vit->property().min_cost + eit->property().cost;
            eit->property().reduced_cost = reduced_cost;
            u_vit->property().sorted_edges_of_vertex.push_back(pair_IntFlt(v,reduced_cost));     // pair: id (outedge), reduced_cost
        }
        sort(u_vit->property().sorted_edges_of_vertex.begin(), u_vit->property().sorted_edges_of_vertex.end(),
            [](pair_IntFlt vecInput_1, pair_IntFlt vecInput_2) {return (vecInput_1.second < vecInput_2.second);} );
#ifdef TOPKSP_DEBUG
        cout<<"the sorted_edges_of_vertex based on reduced cost is:"<<endl;
        cout<<"current ID is "<<u_vit->id()<<". The sorted_edges_of_vertex is: ";
        displayVectorOrListOfPairContents(u_vit->property().sorted_edges_of_vertex); //debug: passed 
#endif        
    }

   
    //// (4) construct the pseudo-tree T
    graph_tau tau;  //T.
    priority_queue<pair_FltInt, vector<pair_FltInt>, min_comp_FltInt> PQ_KSP_candidates_tau; // X.  only store the minCost and interenl ID of tau (the last id). 
    vector<size_t> KSPaths_lastID_tau; // for output, store the top k shortest path id of tau.


    // add the internel_src as the first node in tau
    vertex_iterator_tau src_vit_tau = tau.add_vertex();
    src_vit_tau->property().at_KSPaths = false;//src_vit_tau->property().KSPaths_record.push_back(1); // this node is at the 1st shortest path.
    src_vit_tau->property().internel_id_of_g = internel_src; // this is internel_src
    src_vit_tau->property().predecessor = MY_INFINITY;
    src_vit_tau->property().min_cost = 0; 
    src_vit_tau->property().sum_distance = 0; 
    src_vit_tau->property().sum_hops = 0;
    // construct the first shortest path constructed in tau
    uint64_t internel_src_tau = src_vit_tau->id();
    vertex_iterator_tau dest_vit_tau = add_partialSP_totau(g, internel_src, internel_dest, tau, internel_src_tau);

    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().min_cost, dest_vit_tau->id()) );
    KSPaths_lastID_tau.push_back(dest_vit_tau->id());

 
#ifdef TOPKSP_DEBUG
    cout<<endl;
    cout<<"the inital path (reversed view) at k=1 is "; 
    cout<<"(min_cost,sum_distance,sum_hops)=("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";
    int tmpId = dest_vit_tau->id();    
    vertex_iterator_tau vit_tau_tmp1 = tau.find_vertex(tmpId); 
    cout<<"-->"<<vit_tau_tmp1->property().internel_id_of_g;
    while (tmpId != internel_src_tau)
    {
        tmpId = vit_tau_tmp1->property().predecessor;  
        vit_tau_tmp1 = tau.find_vertex(tmpId);
        cout<<"-->"<<vit_tau_tmp1->property().internel_id_of_g;
    }
    cout<<endl;
    stack<size_t> pk_vkt_nodes_tmp;
#endif
 
    size_t k = 1;
    size_t pk_last_id, tmp_id; // in tau
    size_t tmp_id_g, cur_deviation_id_g; // in tau    
    while (k<Kvalue && !PQ_KSP_candidates_tau.empty())
    {
        // X = X - {p_k}
        PQ_KSP_candidates_tau.pop();  

        // find the deviation path pk_vkt_nodes, the top of pk_vkt_nodes is the deviation node
        stack<size_t> pk_vkt_nodes;
        tmp_id = KSPaths_lastID_tau.back();
        dest_vit_tau = tau.find_vertex(tmp_id);
        dest_vit_tau->property().at_KSPaths = true; //dest_vit_tau->property().KSPaths_record.push_back(k); 
        while (tmp_id != internel_src_tau) //internel_src_tau
        {
            tmp_id = dest_vit_tau->property().predecessor;

            pk_vkt_nodes.push(tmp_id);
            dest_vit_tau = tau.find_vertex(tmp_id);
            if (dest_vit_tau->property().at_KSPaths == true) break;
/*            {
                size_t deviation_id = tmp_id; // not used for this deviation_id. It is the top() of pk_vkt_nodes
                break; // great here
            }*/
            dest_vit_tau->property().at_KSPaths = true; //dest_vit_tau->property().KSPaths_record.push_back(k); 
        }

#ifdef TOPKSP_DEBUG
        pk_vkt_nodes_tmp = pk_vkt_nodes;
        cout<<"the values of pk_vkt_nodes at k="<<k<<" is ";
        while (!pk_vkt_nodes_tmp.empty())
        {
            cout<<pk_vkt_nodes_tmp.top()<<"("<<tau.find_vertex(pk_vkt_nodes_tmp.top())->property().internel_id_of_g<<")"<<",";          
            pk_vkt_nodes_tmp.pop();
        }
        cout<<endl<<endl;
#endif
        // for each node in deviation path, try to find a candidate
        while (!pk_vkt_nodes.empty())  // for each node at pk_vkt
        {
            // find A(v), which is stored in neighborNodes_g (the first of each pair in the vector) and sorted based on reduced cost  
            tmp_id = pk_vkt_nodes.top();
            pk_vkt_nodes.pop();
            vertex_iterator_tau vit_tau_tmp = tau.find_vertex(tmp_id); 
            cur_deviation_id_g = vit_tau_tmp->property().internel_id_of_g;  // cur_deviation_id_g is the v in MPS
            vertex_iterator tmp_vit = g.find_vertex(cur_deviation_id_g);        
            vector<pair_IntFlt> neighborNodes_g = tmp_vit->property().sorted_edges_of_vertex;
            // g.find_vertex_out_neighborNodes(cur_deviation_id_g, neighborNodes_g);  (this func find_vertex_out_neighborNodes is not used now)

            // find A_Tk_(v)
            vector<uint64_t> neighborNodes_tau, neighborNodes_tau_at_g; // 
            tau.find_vertex_inout_neighborNodes(tmp_id, neighborNodes_tau); // here need both in and out neigbor since the tau is directed graph
            for (size_t idx=0; idx<neighborNodes_tau.size(); ++idx)
            {
                vertex_iterator_tau vit_tau_tmptmp = tau.find_vertex(neighborNodes_tau[idx]);
                neighborNodes_tau_at_g.push_back( vit_tau_tmptmp->property().internel_id_of_g );
            }

#ifdef TOPKSP_DEBUG
            cout<< "current k is "<< k <<", current deviation id_g is "<< cur_deviation_id_g << endl;
            cout<< "neighborNodes_g is ";
            displayVectorOrListOfPairContents(neighborNodes_g);
            cout<< "neighborNodes_tau_at_g is ";
            displayVectorOrListContents(neighborNodes_tau_at_g);
#endif

            // check the first arc in A(v) - A_Tk_(v) at Graph g:   neighborNodes_g - neighborNodes_tau_at_g             
            for (size_t idx=0; idx<neighborNodes_g.size(); ++idx)
            {
                if ( !(find(neighborNodes_tau_at_g.begin(), neighborNodes_tau_at_g.end(), neighborNodes_g[idx].first) != neighborNodes_tau_at_g.end()) )
                {
                    // v and x is the cur_deviation_id_g and neighborNodes_g[idx].first, respectively, in MPS alg.
                    // now add point x and edge v->x in tau  vertex_iterator
#ifdef TOPKSP_DEBUG
                    cout<<"the (v,x) is (" << cur_deviation_id_g <<","<<neighborNodes_g[idx].first<<")"<<endl;
                    //cout<<endl;
#endif

                    vertex_iterator_tau dest_vit_tau = tau.add_vertex();
                    

                    edge_iterator eit_g;
                    bool find_result = g.find_out_edge_2id(cur_deviation_id_g, neighborNodes_g[idx].first, eit_g); // for eit_g->property().cost and eit_g->property().phy_dist
                    assert(find_result);

                    dest_vit_tau->property().at_KSPaths = false;//dest_vit_tau->property().KSPaths_record.push_back(1); // this node is at the 1st shortest path.
                    dest_vit_tau->property().predecessor = tmp_id; // tmp_id = pk_vkt_nodes.top();
                    dest_vit_tau->property().internel_id_of_g = neighborNodes_g[idx].first;   
                    dest_vit_tau->property().min_cost     = vit_tau_tmp->property().min_cost + eit_g->property().cost; //   
                    dest_vit_tau->property().sum_distance = vit_tau_tmp->property().sum_distance + eit_g->property().phy_dist; //   
                    dest_vit_tau->property().sum_hops     = vit_tau_tmp->property().sum_hops + 1;

                    //cout<<"the new added node (idx_g,idx_tau) in Tau is ("<<dest_vit_tau->property().internel_id_of_g<<","<<dest_vit_tau->id()<<")"<<endl;
                    //cout<<"the first 3 (min_cost,sum_distance,sum_hops)=("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";

                    uint64_t tau_tmp_id = dest_vit_tau->id();  // for next iteration use
                    edge_iterator_tau eit_tau;
                    tau.add_edge(tmp_id,tau_tmp_id,eit_tau); // note: put all info at vertex, see dest_vit_tau

                    // from neighborNodes_g[idx].first to the dest in the sssp 
                    dest_vit_tau = add_partialSP_totau(g, neighborNodes_g[idx].first, internel_dest, tau, tau_tmp_id);
                    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().min_cost, dest_vit_tau->id()) );

#ifdef TOPKSP_DEBUG
                    cout<<"the new added path (reversed view) at k="<<k<<" is ";  
cout<<"(min_cost,sum_distance,sum_hops)=("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";
                    tmpId = dest_vit_tau->id();    
                    vit_tau_tmp1 = tau.find_vertex(tmpId); 
                    //cout<<"(tmpId,internel_src_tau) is (" << tmpId<<","<<internel_src_tau<<")"<<endl;
                    cout<<"-->"<<vit_tau_tmp1->property().internel_id_of_g;
                    
                    while (tmpId != internel_src_tau)
                    {
                        tmpId = vit_tau_tmp1->property().predecessor;  
                        vit_tau_tmp1 = tau.find_vertex(tmpId);
                        cout<<"-->"<<vit_tau_tmp1->property().internel_id_of_g;
                    }
                    cout<<endl<<endl;
#endif

                    break;
                } 
            }
        }
        
        k++;
        // find the kth shorest path, add to KSPaths_lastID_tau, and add to KSPaths_record
        pk_last_id = PQ_KSP_candidates_tau.top().second; 
        KSPaths_lastID_tau.push_back(pk_last_id);


#ifdef TOPKSP_DEBUG
        cout<<endl;
        cout<<"k just increased now k="<<k<<", and the pk_last_id="<<pk_last_id<<endl;
#endif

    }


    // output_top_ksp: 
    cout<<"The top KSP Results:"<<endl;
    for (size_t idx = 0; idx < KSPaths_lastID_tau.size(); ++idx)
    {
        size_t tmp_k = idx + 1;
        uint64_t curr_id_tau = KSPaths_lastID_tau[idx];

        vertex_iterator_tau vit_tau = tau.find_vertex(curr_id_tau);
        cout<<endl;
        cout<<"The "<<tmp_k<<" result with metrics:   ";
        cout<<"min_cost is "    <<vit_tau->property().min_cost<<"; ";
        cout<<"sum_distance is "<<vit_tau->property().sum_distance<<"; ";
        cout<<"sum_hops is "    <<vit_tau->property().sum_hops<<"."<<endl;

        uint64_t curr_id_g = vit_tau->property().internel_id_of_g; 
        uint64_t curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
        stack<uint64_t> curr_path;
        curr_path.push(curr_exID);
        do
        {
            curr_id_tau = vit_tau->property().predecessor;
            vit_tau = tau.find_vertex(curr_id_tau);
            curr_id_g = vit_tau->property().internel_id_of_g; 
            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            curr_path.push(curr_exID);
        } while (curr_id_tau != internel_src_tau);
        cout<<"The "<<tmp_k<<" result with path nodes: ";
        while (!curr_path.empty())
        {
            cout<<"-->"<<curr_path.top();
            curr_path.pop();
        }
        cout<<endl;
    }
 

#ifdef SIM
    SIM_END(true);
#endif
    perf.stop(perf_group);
    return;
}


 


/*void sssp(graph_t& g, size_t src, gBenchPerf_event & perf, int perf_group)  // src is external ID
{
    priority_queue<pair_IntFlt, vector<pair_IntFlt>, comp> PQ;
    
    perf.open(perf_group);
    perf.start(perf_group);
#ifdef SIM
    SIM_BEGIN(true);
#endif
    // initialize
    uint64_t internel_src = g.external_to_internel_id(to_string(src)); // src is external id, all the processing is related internel id.

    vertex_iterator src_vit = g.find_vertex(internel_src); // src
    src_vit->property().min_cost = 0;
    PQ.push(pair_IntFlt(internel_src,0)); // src

    // for every un-visited vertex, try relaxing the path
    while (!PQ.empty())
    {
        size_t u = PQ.top().first; 
        PQ.pop();

        vertex_iterator u_vit = g.find_vertex(u);

        for (edge_iterator eit = u_vit->edges_begin(); eit != u_vit->edges_end(); eit++)
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);

            //size_t alt = u_vit->property().min_cost + eit->property().cost; //cost  phy_dist
            float alt = u_vit->property().min_cost + eit->property().cost; //cost  phy_dist
            if (alt < v_vit->property().min_cost) 
            {
                v_vit->property().min_cost = alt;
                v_vit->property().predecessor = u;
                PQ.push(pair_IntFlt(v,alt));
            }
        }
    }
#ifdef SIM
    SIM_END(true);
#endif
    perf.stop(perf_group);
    return;
}*/
/*void print_path(const graph_t& g, size_t src, size_t dest) // src, dest, and processing are based on internel ID. print out is the external ID.
{
    vertex_iterator vit = g.find_vertex(dest);
    uint64_t dest_parent = vit->property().predecessor;
    if (dest == src)
        cout<<"-->"<< g.internal_to_externel_id(src);
    else if (dest_parent == MY_INFINITY)
        cout<<"no path from "<< g.internal_to_externel_id(src) << " to " << g.internal_to_externel_id(dest) <<" exists"<< endl;
    else
        print_path(g, src, dest_parent);
    cout<<"-->"<< g.internal_to_externel_id(dest);
}
*/
/*void print_path(graph_t& g, size_t src, size_t dest) // src, dest are external ID. // and processing are based on internel ID. print out is the external ID.
{
    uint64_t internel_dest = g.external_to_internel_id(to_string(dest));
    vertex_iterator vit = g.find_vertex(internel_dest);
    uint64_t dest_parent = vit->property().predecessor;
    if (dest == src)
        cout<<"-->"<< src;
    else if (dest_parent == MY_INFINITY)
        cout<<"no path from "<< src << " to " << dest <<" exists"<< endl;
    else
    {
        print_path(g, src, atoi(g.internal_to_externel_id(dest_parent).c_str()));
        cout<<"-->"<< dest;
    }
    
}
//==============================================================//
void output(graph_t& g, uint64_t root) // root is exID
{
    cout<<"Results: \n";
    vertex_iterator vit;
    for (vit=g.vertices_begin(); vit!=g.vertices_end(); vit++) // obtained is intenel ID
    {
        //cout<<"== vertex "<<vit->id()<<": min_cost-";    //std::string internal_to_externel_id(uint64_t inID) // new
        uint64_t curr_exID = atoi(g.internal_to_externel_id(vit->id()).c_str());
        cout<<"vertex "<<root<<" to vertex "<< curr_exID <<": MinCost is ";  // output the corresponding externel id
        if (vit->property().min_cost == MY_INFINITY) 
            cout<<"INF";
        else
        {
            cout<<vit->property().min_cost;
            cout<<". The path by exID is: ";
            print_path(g, root, curr_exID);
        }

        cout<<"\n";
    }
    return;
}*/






void reset_graph(graph_t & g)
{
    vertex_iterator vit;
    for (vit=g.vertices_begin(); vit!=g.vertices_end(); vit++)
    {
        vit->property().successor = MY_INFINITY;
        vit->property().min_cost = FLT_MAX;
        vit->property().sorted_edges_of_vertex.clear();
    }

}
//==============================================================//
//==============================================================//
void arg_init(argument_parser & arg)
{
    //arg.add_arg("root","0","root/starting vertex");
    arg.add_arg("src","0","root/src vertex");
    arg.add_arg("dest","0","root/dest vertex");
    arg.add_arg("Kvalue","1","root/K value");
}
//==============================================================//






int main(int argc, char * argv[])
{
    graphBIG::print();
    cout<<"Benchmark: sssp shortest path\n";
    double t1, t2;

    argument_parser arg;
    gBenchPerf_event perf;
    arg_init(arg);
    if (arg.parse(argc,argv,perf,false)==false)
    {
        arg.help();
        return -1;
    }
    string path, separator;
    arg.get_value("dataset",path);
    arg.get_value("separator",separator);

/*    size_t root,threadnum;
    arg.get_value("root",root);
    arg.get_value("threadnum",threadnum);*/

    size_t src, dest, Kvalue, threadnum;
    arg.get_value("src",src);
    arg.get_value("dest",dest);
    arg.get_value("Kvalue",Kvalue);
    arg.get_value("threadnum",threadnum);



    graph_t graph;
    cout<<"loading data... \n";

    t1 = timer::get_usec();
    string vfile = path + "/vertex.csv";
    string efile = path + "/edge.csv";

#ifndef EDGES_ONLY    
    if (graph.load_csv_vertices(vfile, true, separator, 0) == -1)
        return -1;
    if (graph.load_csv_edges(efile, true, separator, 0, 1) == -1) 
        return -1;
#else
    if (graph.load_csv_edges(efile, true, separator, 0, 1) == -1)
        return -1;
#endif

    size_t vertex_num = graph.num_vertices();
    size_t edge_num = graph.num_edges();
    t2 = timer::get_usec();
    cout<<"== "<<vertex_num<<" vertices  "<<edge_num<<" edges\n";
#ifndef ENABLE_VERIFY
    cout<<"== time: "<<t2-t1<<" sec\n\n";
#endif


    // sanity check
    uint64_t internel_srcID = graph.external_to_internel_id(to_string(src));
    if (graph.find_vertex(internel_srcID)==graph.vertices_end()) 
 
    {
        cerr<<"wrong source vertex: "<<src<<endl;
        return 0;
    }
    uint64_t internel_destID = graph.external_to_internel_id(to_string(dest));
    if (graph.find_vertex(internel_destID)==graph.vertices_end()) 
    {
        cerr<<"wrong dest vertex: "<<dest<<endl;
        return 0;
    }

 
    //cout<<"Shortest Path: source-"<<root;
    cout<<"...\n";

    gBenchPerf_multi perf_multi(threadnum, perf);
    unsigned run_num = ceil(perf.get_event_cnt() /(double) DEFAULT_PERF_GRP_SZ);
    if (run_num==0) run_num = 1;
    double elapse_time = 0;
    


    //cout<<"src is "<<src<<" dest is "<<dest<<" Kvalue is "<<Kvalue<<endl;
    //std::cout<<"the run_num is="<<run_num==1<<std::endl; // tested, run_num==1
    for (unsigned i=0;i<run_num;i++)
    {
        t1 = timer::get_usec();

        if (threadnum==1)
            //sssp(graph, root, perf, i);
            top_ksp(graph, src, dest, Kvalue, perf, i);
        else
            //parallel_sssp(graph, root, threadnum, perf_multi, i);
        
        t2 = timer::get_usec();
        elapse_time += t2-t1;
        if ((i+1)<run_num) reset_graph(graph);
    }
#ifndef ENABLE_VERIFY
    cout<<"== time: "<<elapse_time/run_num<<" sec\n";
    if (threadnum == 1)
        perf.print();
    else
        perf_multi.print();
#endif


/*#ifdef ENABLE_OUTPUT
    cout<<"\n";
    output(graph,root);
#endif*/

    cout<<"==================================================================\n";
    return 0;
}  // end main

