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
#include "omp.h"
#include <queue>
#include <cfloat>     
#include <cassert>
#include <stack>
#include <cstdlib>
using namespace std;

//#define TOPKSP_PRINTOUT
//#define VERIFY_RESULTS
#define MY_INFINITY  0xfff0  // for unsigned i.e., 65520
#define ENABLE_OUTPUT       //new

#ifdef HMC
#include "HMC.h"
#endif
#ifdef SIM
#include "SIM.h"
#endif

size_t beginiter = 0;
size_t enditer = 0;

typedef pair<size_t,float> pair_IntFlt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<float,size_t> pair_FltInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //


class vertex_property   // new
{
public:
    vertex_property():min_cost(FLT_MAX),successor(MY_INFINITY){}
    float min_cost;     // for shortest path
    uint64_t successor; // for shortest path 
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
        vit->property().successor = MY_INFINITY;
        vit->property().min_cost = FLT_MAX;
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
    vertex_property_tau():at_KSPaths(false),predecessor(MY_INFINITY),internal_id_of_g(MY_INFINITY),min_cost(0),sum_distance(0),sum_hops(0){}

    bool at_KSPaths;           //vector<size_t> KSPaths_record;
    uint64_t predecessor;      // for shortest path, store the id of the graph tau
    uint64_t internal_id_of_g; // internal id of orginal graph g

    //for shortest path
    float min_cost;            
    float sum_distance;        
    size_t sum_hops;           
};


class edge_property_tau // new
{
public:
    edge_property_tau(){} //edge_property_tau():cost(FLT_MAX),phy_dist(FLT_MAX),reduced_cost(FLT_MAX){} 
};

typedef openG::extGraph<vertex_property_tau, edge_property_tau> graph_tau;
typedef graph_tau::vertex_iterator    vertex_iterator_tau;
typedef graph_tau::edge_iterator      edge_iterator_tau;

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
        uint64_t internal_dest = g.external_to_internel_id(to_string(dest)); 
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
        } while (curr_id_g != internal_dest);
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
    assert(src_vit_tau->property().internal_id_of_g == src_id_g); // this is needed for the following code to work

    uint64_t tau_tmp_id = start_id_tau;
    src_vit_g = g.find_vertex(src_id_g);  
    while (src_vit_g->id() != dest_id_g) // note: dest_id_g should be the internal_dest
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
        dest_vit_tau->property().internal_id_of_g = dest_vit_g->id();   
        dest_vit_tau->property().min_cost     = src_vit_tau->property().min_cost + eit_g->property().cost; //   
        dest_vit_tau->property().sum_distance = src_vit_tau->property().sum_distance + eit_g->property().phy_dist; //   
        dest_vit_tau->property().sum_hops     = src_vit_tau->property().sum_hops + 1;

        //cout<<"the new added node (idx_g,idx_tau) in Tau is ("<<dest_vit_tau->property().internal_id_of_g<<","<<dest_vit_tau->id()<<")"<<endl;
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

bool is_loopless_path(graph_t& g, graph_tau& tau, size_t path_last_id_tau, size_t path_first_id_tau)//based on sorted vector
{
    size_t tmpId = path_last_id_tau;    
    vertex_iterator_tau vit_tau_tmp = tau.find_vertex(tmpId); 

    vector<size_t> path_id_set_g;  
    path_id_set_g.push_back( vit_tau_tmp->property().internal_id_of_g );
    while (tmpId != path_first_id_tau)
    {
        tmpId = vit_tau_tmp->property().predecessor;  
        vit_tau_tmp = tau.find_vertex(tmpId);
        path_id_set_g.push_back( vit_tau_tmp->property().internal_id_of_g );
    }
    sort(path_id_set_g.begin(), path_id_set_g.end());
    return adjacent_find(path_id_set_g.begin(), path_id_set_g.end()) == path_id_set_g.end();
}


void top_ksp(graph_t& g, size_t src, size_t dest, size_t Kvalue, double max_phy_dist, size_t max_phy_hops, gBenchPerf_event & perf, int perf_group, ofstream &res_fstream) // src and dest are exID
{
    perf.open(perf_group);
    perf.start(perf_group);
#ifdef SIM
    SIM_BEGIN(true);
#endif
    uint64_t internal_src  = g.external_to_internel_id(to_string(src)); 
    uint64_t internal_dest = g.external_to_internel_id(to_string(dest)); 

    //// Now all processing is based on internal id first.
    //// (1) sssp sub-procedure  -->for vertex: update  v_vit->property().min_cost and v_vit->property().predecessor 
    priority_queue<pair_IntFlt, vector<pair_IntFlt>, min_comp_IntFlt> PQ;

    vertex_iterator dest_vit = g.find_vertex(internal_dest);  // note: here the source of sssp is internal_dest, 
    dest_vit->property().min_cost = 0;
    PQ.push(pair_IntFlt(internal_dest,0));

    // vit->property().predecessor is used to construct sssp, where the ancestor of all nodes is internal_dest, 
    // by using vit->property().predecessor recursively, we can find the shortest path of any node to internal_dest
    while (!PQ.empty()) 
    {
        size_t u = PQ.top().first; //id
        PQ.pop();

        vertex_iterator u_vit = g.find_vertex(u);
        for (edge_iterator eit = u_vit->in_edges_begin(); eit != u_vit->in_edges_end(); eit++)  // new for in edge
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);
            // for every vertex, try relaxing the path
            float alt = u_vit->property().min_cost + eit->property().cost; // 
            if (alt < v_vit->property().min_cost) 
            {
                v_vit->property().min_cost = alt; 
                v_vit->property().successor = u;    // new, ori is predecessor
                PQ.push(pair_IntFlt(v,alt));
            }
        }
    }
    
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

    // add the internal_src as the first node in tau
    vertex_iterator_tau src_vit_tau = tau.add_vertex();
    src_vit_tau->property().at_KSPaths = false;//src_vit_tau->property().KSPaths_record.push_back(1); // this node is at the 1st shortest path.
    src_vit_tau->property().internal_id_of_g = internal_src; // this is internal_src
    src_vit_tau->property().predecessor = MY_INFINITY;
    src_vit_tau->property().min_cost = 0; 
    src_vit_tau->property().sum_distance = 0; 
    src_vit_tau->property().sum_hops = 0;
    // construct the first shortest path constructed in tau
    uint64_t internal_src_tau = src_vit_tau->id();

#ifdef TOPKSP_DEBUG
    cout<< "I come before add_partialSP_totau" <<endl;
    cout << "internal_src = "<< internal_src <<endl;
    cout << "internal_dest = " << internal_dest <<endl;
    cout << "internal_src_tau = " << internal_src_tau << endl;
    
    // output externel ID
    for (int src_tmp=0; src_tmp < g.num_vertices()-1; src_tmp++)//for MPS alg simple ex. note: cannot output dest to dest
    {
        if (src_tmp != dest)
        {
            cout<<"the shortest path: "<<src_tmp<<"----->"<<dest<<":"<<endl;
            output_shorest_path(g, src_tmp, dest);
        }

    }//may terminate if there is no path between to nodes
#endif  

    vertex_iterator_tau dest_vit_tau = add_partialSP_totau(g, internal_src, internal_dest, tau, internal_src_tau);
    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().min_cost, dest_vit_tau->id()) );


 
#ifdef TOPKSP_DEBUG
    cout<<endl;
    cout<<"the inital path (reversed view) at k=1 is "; 
    cout<<"(min_cost,sum_distance,sum_hops)=("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";
    int tmpId = dest_vit_tau->id();    
    vertex_iterator_tau vit_tau_tmp1 = tau.find_vertex(tmpId); 
    cout<<"-->"<<vit_tau_tmp1->property().internal_id_of_g;
    while (tmpId != internal_src_tau)
    {
        tmpId = vit_tau_tmp1->property().predecessor;  
        vit_tau_tmp1 = tau.find_vertex(tmpId);
        cout<<"-->"<<vit_tau_tmp1->property().internal_id_of_g;
    }
    cout<<endl;
    stack<size_t> pk_vkt_nodes_tmp;
#endif
 

    size_t k = 0;
    size_t iter = 0; 
    const size_t max_iter = Kvalue*10;
    while (k<Kvalue && !PQ_KSP_candidates_tau.empty() && iter<max_iter)
    {
        // find the current shortest path p (may have loop)
        size_t candi_last_id = PQ_KSP_candidates_tau.top().second; 

        // X = X - {p}
        PQ_KSP_candidates_tau.pop();          

        // if p is loopless and within max_phy_dist and max_phy_hops, add it to the final KSP,   
        bool is_loopless = is_loopless_path(g, tau, candi_last_id, internal_src_tau);
        bool within_max_phy_dist = tau.find_vertex(candi_last_id)->property().sum_distance <= max_phy_dist;
        bool within_max_phy_hops = tau.find_vertex(candi_last_id)->property().sum_hops     <= max_phy_hops;

        if (is_loopless && within_max_phy_dist && within_max_phy_hops)  
        {
            KSPaths_lastID_tau.push_back(candi_last_id);
            k++;
        }

        //note: even the current shortest path p may have loop, the new generated candidates based on p may have no loop.
        // find the deviation path pk_vkt_nodes, the top of pk_vkt_nodes is the deviation node
        stack<size_t> pk_vkt_nodes;
        size_t tmp_id = candi_last_id; // KSPaths_lastID_tau.back()
        dest_vit_tau = tau.find_vertex(tmp_id);
        dest_vit_tau->property().at_KSPaths = true;  
        while (tmp_id != internal_src_tau)  
        {
            tmp_id = dest_vit_tau->property().predecessor;

            pk_vkt_nodes.push(tmp_id);
            dest_vit_tau = tau.find_vertex(tmp_id);
            if (dest_vit_tau->property().at_KSPaths == true) break;
/*            {
                size_t deviation_id = tmp_id; // not used for this deviation_id. It is the top() of pk_vkt_nodes
                break; // great here
            }*/
            dest_vit_tau->property().at_KSPaths = true;  
        }

#ifdef TOPKSP_DEBUG
        pk_vkt_nodes_tmp = pk_vkt_nodes;
        cout<<"the values of pk_vkt_nodes at iter="<<iter<<" is ";
        while (!pk_vkt_nodes_tmp.empty())
        {
            cout<<pk_vkt_nodes_tmp.top()<<"("<<tau.find_vertex(pk_vkt_nodes_tmp.top())->property().internal_id_of_g<<")"<<",";          
            pk_vkt_nodes_tmp.pop();
        }
        cout<<endl<<endl;
#endif
        // for each node in deviation path, try to find a candidate
        while (!pk_vkt_nodes.empty())  // for each node at pk_vkt
        {
            // find current deviation point: pk_top_id
            size_t pk_top_id = pk_vkt_nodes.top();  // tmp_id--> pk_top_id
            pk_vkt_nodes.pop();

            // get out if there is loop to save time
            if ( !is_loopless_path(g, tau, pk_top_id, internal_src_tau) )
            {
                break;
            }

            // find A(v), which is stored in neighborNodes_g (the first of each pair in the vector) and sorted based on reduced cost  
            vertex_iterator_tau pk_top_vit = tau.find_vertex(pk_top_id);   //vit_tau_tmp--> pk_top_vit
            size_t cur_deviation_id_g = pk_top_vit->property().internal_id_of_g;  // cur_deviation_id_g is the v in MPS
            vertex_iterator tmp_vit = g.find_vertex(cur_deviation_id_g);        
            vector<pair_IntFlt> neighborNodes_g = tmp_vit->property().sorted_edges_of_vertex;
            // g.find_vertex_out_neighborNodes(cur_deviation_id_g, neighborNodes_g);  (this func find_vertex_out_neighborNodes is not used now)

/*            // find A_Tk_(v),  (the old version for the case with loop)
            vector<uint64_t> neighborNodes_tau, neighborNodes_tau_at_g; // 
            tau.find_vertex_inout_neighborNodes(pk_top_id, neighborNodes_tau); // here need both in and out neigbor since the tau is directed graph
            for (size_t idx=0; idx<neighborNodes_tau.size(); ++idx)
            {
                vertex_iterator_tau vit_tau_tmptmp = tau.find_vertex(neighborNodes_tau[idx]);
                neighborNodes_tau_at_g.push_back( vit_tau_tmptmp->property().internal_id_of_g );
            }*/

            // (new version )find A_Tk_(v),  now further update it that  A_Tk_(v) contains out edges (sucsessors) and all predecessors up to internal_src_tau
            // this update is to reduce the number of loop paths
            vector<uint64_t> neighborNodes_tau, neighborNodes_tau_at_g; // 
            tau.find_vertex_out_neighborNodes(pk_top_id, neighborNodes_tau); // here only need  out neigbor, the
            for (size_t idx=0; idx<neighborNodes_tau.size(); ++idx)
            {
                vertex_iterator_tau vit_tau_tmptmp = tau.find_vertex(neighborNodes_tau[idx]);
                neighborNodes_tau_at_g.push_back( vit_tau_tmptmp->property().internal_id_of_g );
            }
            tmp_id = pk_top_id;
            vertex_iterator_tau vit_tau_tmp = pk_top_vit;
            while (tmp_id != internal_src_tau)
            {
                tmp_id = vit_tau_tmp->property().predecessor;  
                vit_tau_tmp = tau.find_vertex(tmp_id);
                neighborNodes_tau_at_g.push_back(vit_tau_tmp->property().internal_id_of_g);
            }


#ifdef TOPKSP_DEBUG
            cout<< "current iter is "<< iter <<", current deviation id_g is "<< cur_deviation_id_g << endl;
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
                    dest_vit_tau->property().predecessor = pk_top_id; // pk_top_id = pk_vkt_nodes.top();
                    dest_vit_tau->property().internal_id_of_g = neighborNodes_g[idx].first;   
                    dest_vit_tau->property().min_cost     = pk_top_vit->property().min_cost + eit_g->property().cost; //   
                    dest_vit_tau->property().sum_distance = pk_top_vit->property().sum_distance + eit_g->property().phy_dist; //   
                    dest_vit_tau->property().sum_hops     = pk_top_vit->property().sum_hops + 1;

                    //cout<<"the new added node (idx_g,idx_tau) in Tau is ("<<dest_vit_tau->property().internal_id_of_g<<","<<dest_vit_tau->id()<<")"<<endl;
                    //cout<<"the first 3 (min_cost,sum_distance,sum_hops)=("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";

                    uint64_t tau_tmp_id = dest_vit_tau->id();  // for next iteration use
                    edge_iterator_tau eit_tau;
                    tau.add_edge(pk_top_id,tau_tmp_id,eit_tau); // note: put all info at vertex, see dest_vit_tau

                    // from neighborNodes_g[idx].first to the dest in the sssp 


                    //cout << "neighborNodes_g[idx].first = "<< neighborNodes_g[idx].first<<endl;
                    //cout << "internal_dest = " << internal_dest <<endl;
                    //cout << "tau_tmp_id = " << tau_tmp_id << endl;

                    if (neighborNodes_g[idx].first != internal_dest)
                    {
                        dest_vit_tau = add_partialSP_totau(g, neighborNodes_g[idx].first, internal_dest, tau, tau_tmp_id);
                    }    
                    

                    //cout<< "I come here before PQ_KSP_candidates_tau"<<endl;
                    //cout<<"dest_vit_tau->property().min_cost" << dest_vit_tau->property().min_cost << endl;
                    //cout<<"dest_vit_tau->id()" << dest_vit_tau->id() << endl;

                    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().min_cost, dest_vit_tau->id()) );


#ifdef TOPKSP_DEBUG
                    cout<<"the new added path (reversed view) at iter="<<iter<<" is ";  
cout<<"(min_cost,sum_distance,sum_hops)=("<<dest_vit_tau->property().min_cost<<","<<dest_vit_tau->property().sum_distance<<","<<dest_vit_tau->property().sum_hops<<")";
                    tmpId = dest_vit_tau->id();    
                    vit_tau_tmp1 = tau.find_vertex(tmpId); 
                    //cout<<"(tmpId,internal_src_tau) is (" << tmpId<<","<<internal_src_tau<<")"<<endl;
                    cout<<"-->"<<vit_tau_tmp1->property().internal_id_of_g;
                    
                    while (tmpId != internal_src_tau)
                    {
                        tmpId = vit_tau_tmp1->property().predecessor;  
                        vit_tau_tmp1 = tau.find_vertex(tmpId);
                        cout<<"-->"<<vit_tau_tmp1->property().internal_id_of_g;
                    }
                    cout<<endl<<endl;
#endif

                    break;
                } 
            } 
        }//while (!pk_vkt_nodes.empty())  // for each node at pk_vkt
        
        iter++;



#ifdef TOPKSP_DEBUG
        cout<<endl;
        cout<<"iter just increased now iter="<<iter<<endl;//", and the pk_last_id="<<pk_last_id<<endl;
#endif

    } //while (k<Kvalue && !PQ_KSP_candidates_tau.empty() && iter<max_iter)



#ifdef TOPKSP_PRINTOUT
    // output_top_ksp: 
    res_fstream<<"The top "<<Kvalue<<" shortest loopless paths";
    res_fstream<<"(the acutal number of running iteration is "<<iter<<"):"<<endl;
    if (KSPaths_lastID_tau.size()==0)
    {
        res_fstream<<"cannot find any path, maybe we have too much constraints!"<<endl;
        assert(false);
    }
    for (size_t idx = 0; idx < KSPaths_lastID_tau.size(); ++idx)
    {
        size_t tmp_k = idx + 1;
        uint64_t curr_id_tau = KSPaths_lastID_tau[idx];

        vertex_iterator_tau vit_tau = tau.find_vertex(curr_id_tau);
        res_fstream<<endl;
        res_fstream<<"The "<<tmp_k<<" result with metrics:   ";
        res_fstream<<"min_cost is "    <<vit_tau->property().min_cost<<"; ";
        res_fstream<<"sum_distance is "<<vit_tau->property().sum_distance<<"; ";
        res_fstream<<"sum_hops is "    <<vit_tau->property().sum_hops<<"."<<endl;

        uint64_t curr_id_g = vit_tau->property().internal_id_of_g; 
        uint64_t curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
        stack<uint64_t> curr_path;
        curr_path.push(curr_exID);
        do
        {
            curr_id_tau = vit_tau->property().predecessor;
            vit_tau = tau.find_vertex(curr_id_tau);
            curr_id_g = vit_tau->property().internal_id_of_g; 
            curr_exID = atoi(g.internal_to_externel_id(curr_id_g).c_str());
            curr_path.push(curr_exID);
        } while (curr_id_tau != internal_src_tau);
        res_fstream<<"The "<<tmp_k<<" result with path nodes: ";
        while (!curr_path.empty())
        {
            res_fstream<<"-->"<<curr_path.top();
            curr_path.pop();
        }

        if ( !is_loopless_path(g, tau, KSPaths_lastID_tau[idx], internal_src_tau) )
        {
            res_fstream<<endl<<"wrong, this path has loop"<<endl;
            assert(false);
        }
        res_fstream<<endl;
    }
    if (Kvalue!=KSPaths_lastID_tau.size())
    {
        res_fstream<<"Warning: the input k= "<<Kvalue<<" is too large, cannot find the required loopless paths with current constraints!"<<endl;
        res_fstream<<"the largest possible k is "<<KSPaths_lastID_tau.size()<<endl;
    }
#endif


#ifdef SIM
    SIM_END(true);
#endif
    perf.stop(perf_group);
    return;
}


//==============================================================//
//==============================================================//
void arg_init(argument_parser & arg)
{
    //arg.add_arg("root","0","root/starting vertex");
    arg.add_arg("src","0","root/src vertex");
    arg.add_arg("dest","10","root/dest vertex");
    arg.add_arg("Kvalue","3","root/K value");
    arg.add_arg("max_phy_dist","100000","root/max_phy_dist");
    arg.add_arg("max_phy_hops","100000","root/max_phy_hops");
}
//==============================================================//


//Initialize a graph
void graphInit(graph_t &graph, string vfile, string efile, size_t src, size_t dest, string separator)
{
    cout<<"loading data... \n";
    double t1 = timer::get_usec();
    
    #ifndef EDGES_ONLY    
        if (graph.load_csv_vertices(vfile, true, separator, 0) == -1)
            return;
        if (graph.load_csv_edges(efile, true, separator, 0, 1) == -1) 
            return;
    #else
        if (graph.load_csv_edges(efile, true, separator, 0, 1) == -1)
            return;
    #endif

    size_t vertex_num = graph.num_vertices();
    size_t edge_num = graph.num_edges();
    double t2 = timer::get_usec();
    cout<<"== "<<vertex_num<<" vertices  "<<edge_num<<" edges\n";

    #ifndef ENABLE_VERIFY
        cout<<"== time: "<<t2-t1<<" sec\n\n";
    #endif

    // sanity check
    uint64_t internel_srcID = graph.external_to_internel_id(to_string(src));
    if (graph.find_vertex(internel_srcID)==graph.vertices_end()) 
 
    {
        cerr<<"wrong source vertex: "<<src<<endl;
        assert(false);
    }
    uint64_t internel_destID = graph.external_to_internel_id(to_string(dest));
    if (graph.find_vertex(internel_destID)==graph.vertices_end()) 
    {
        cerr<<"wrong dest vertex: "<<dest<<endl;
        assert(false);
    }

}

//Generate several test pairs
vector<vector<int> > genTests(int total_vertex, int AGG_pair_num)
{
    vector<vector<int> > res;
    int min, max, p_src, p_dest;
    srand (1223);

    for (int i=0; i<AGG_pair_num; i++)
    {        
        min = 0;
        max = 19;
        p_src = rand()%(max-min + 1) + min;
        
        min  = total_vertex-20;
        max  = total_vertex-1;                
        p_dest = rand()%(max-min + 1) + min;
        vector<int> pair = {p_src, p_dest}; 
        res.push_back(pair);
    }
    return res;
}


double serialTest(string vfile, string efile, size_t src, size_t dest, string separator, int AGG_pair_num, int Kvalue, int max_phy_dist, int max_phy_hops, vector<vector<int> >tests, gBenchPerf_event perf)
{
    graph_t *graph;  
    double t1 = timer::get_usec();
    int p_src, p_dest;

    for (int i=0; i<AGG_pair_num; i++)
    {
        graph = new graph_t();
        graphInit(*graph, vfile, efile, src, dest, separator);
        p_src = tests[i][0]; 
        p_dest = tests[i][1];
        ofstream res_fstream;

        #ifdef TOPKSP_PRINTOUT
            res_fstream.open("./test_res/serial_res_"+to_string(i), ofstream::trunc);
            res_fstream<<endl<<endl<<"Now output the top "<<Kvalue<<" Shortest paths bwteen nodes: "<<p_src<<"-->"<<p_dest;
            res_fstream<<"(run "<<i+1<<" out of "<<AGG_pair_num<<")"<<endl;
        #endif
        
        top_ksp(*graph, p_src , p_dest, Kvalue, max_phy_dist, max_phy_hops, perf, 0, res_fstream);  
        delete graph;

        #ifdef TOPKSP_PRINTOUT
            res_fstream.close();
        #endif
    }

    double t2 = timer::get_usec();
    return t2-t1;
}



double parallelTest(string vfile, string efile, size_t src, size_t dest, string separator, int AGG_pair_num, int Kvalue, int max_phy_dist, int max_phy_hops, vector<vector<int> >tests, gBenchPerf_event perf)
{
    graph_t *graph;     
    double t1 = timer::get_usec();

    omp_set_num_threads(AGG_pair_num);
    #pragma omp parallel private(graph)
    {
        int idx = omp_get_thread_num();
        int p_src = tests[idx][0];
        int p_dest = tests[idx][1];
        ofstream res_fstream;

        graph = new graph_t();
        graphInit(*graph, vfile, efile, src, dest, separator);

        #ifdef TOPKSP_PRINTOUT
            res_fstream.open("./test_res/parallel_res_"+to_string(idx), ofstream::trunc);
            res_fstream<<endl<<endl<<"Now output the top "<<Kvalue<<" Shortest paths bwteen nodes: "<<p_src<<"-->"<<p_dest;
            res_fstream<<"(run "<<idx+1<<" out of "<<AGG_pair_num<<")"<<endl;
        #endif

        top_ksp(*graph, p_src , p_dest, Kvalue, max_phy_dist, max_phy_hops, perf, 0, res_fstream);  
        delete graph;

        #ifdef TOPKSP_PRINTOUT
            res_fstream.close();
        #endif
    }

    double t2 = timer::get_usec(); 
    return t2-t1;
}


int main(int argc, char * argv[])
{
    graphBIG::print();
    cout<<"Benchmark: sssp shortest path\n";

    argument_parser arg;
    gBenchPerf_event perf;
    arg_init(arg);
    if (arg.parse(argc,argv,perf,false)==false)
    {
        arg.help();
        return -1;
    }
    string path, separator;
    size_t src, dest, Kvalue, max_phy_hops, threadnum, AGG_pair_num, total_vertex;
    double max_phy_dist;
    string vfile, efile;

    arg.get_value("dataset",path);
    arg.get_value("separator",separator);
    arg.get_value("src",src);
    arg.get_value("dest",dest);
    arg.get_value("Kvalue",Kvalue);
    arg.get_value("max_phy_dist",max_phy_dist);
    arg.get_value("max_phy_hops",max_phy_hops);
    arg.get_value("threadnum",threadnum);

    
    AGG_pair_num = 24;   // 5    50    100 
    Kvalue = 30;         // 5   50    500  
    max_phy_dist = 10;
    max_phy_hops = 8;
    total_vertex = 400;  // 400   4000   10000 
    vfile = path + "/vertex400.csv";  // vertex400   vertex4000    vertex10000
    efile = path + "/vertex400dim8_edge.csv";       
    //  vertex400dim8_edge     vertex400dim12_edge       vertex400dim16_edge
    //  vertex4000dim8_edge    vertex4000dim12_edge      vertex4000dim16_edge
    //  vertex10000dim8_edge   vertex10000dim12_edge     vertex10000dim16_edge

    //cout<<"Shortest Path: source-"<<root;
    cout<<"...\n";

    gBenchPerf_multi perf_multi(threadnum, perf);

    //Generate tests
    vector<vector<int> > tests = genTests(total_vertex, AGG_pair_num);
    double s_time, p_time;

    cout<<"AGG_pair_num is "<< AGG_pair_num <<"; Kvalue is "<<Kvalue<<endl;
    cout<<"max_phy_dist is " << max_phy_dist <<endl;
    cout<<"max_phy_hops is " << max_phy_hops <<endl<<endl;
    
    //Serial run
    cout<<"Start running serial test"<<endl;
    s_time = serialTest(vfile, efile, src, dest, separator, AGG_pair_num, Kvalue, max_phy_dist, max_phy_hops, tests, perf);

    #ifndef ENABLE_VERIFY
        cout<<"== Total ruuning time: "<<s_time<<" sec\n";
        cout<<"==================================================================\n"<<endl;
    #endif    

    //Parallel run
    cout<<"Start running openMP test"<<endl;
    p_time = parallelTest(vfile, efile, src, dest, separator, AGG_pair_num, Kvalue, max_phy_dist, max_phy_hops, tests, perf);
    
    #ifndef ENABLE_VERIFY
        cout<<"== Total ruuning time: "<<p_time<<" sec\n";
        cout<<"==================================================================\n"<<endl;
    #endif    

    //Reseult verification
    #ifdef VERIFY_RESULTS
        bool all_identical = 1;
        for (unsigned int i=0; i<AGG_pair_num; i++)
        {
            string cmd = "diff ./test_res/serial_res_"+to_string(i)+" ./test_res/parallel_res_"+to_string(i);
            bool is_identical = system(cmd.c_str())!=0;
            all_identical &= (!is_identical);
            if (is_identical)
                cout<<"Results differ for AGG_pair_"<<i<<endl;
        }
        if (all_identical)
            cout<<"All results are identical!"<<endl;
    #endif

    cout<<"==Speed up: "<<s_time/p_time<<endl;
    return 0;
}  // end main

