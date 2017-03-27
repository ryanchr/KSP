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

#ifdef HMC
#include "HMC.h"
#endif

#ifdef SIM
#include "SIM.h"
#endif

#define OUTPUT_PATHS
#define TOPKSP_PRINTOUT
//#define VERIFY_RESULTS
#define MY_INFINITY  0xfff0  // for unsigned i.e., 65520
#define ENABLE_OUTPUT       //new

size_t beginiter = 0;
size_t enditer = 0;

typedef pair<size_t,float> pair_IntFlt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<float,size_t> pair_FltInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<size_t,size_t> pair_IntInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //

class vertex_property // new
{
public:  // sum_distance(0),sum_hops(0){}
    vertex_property():min_cost(FLT_MAX),successor(MY_INFINITY),sum_distance(FLT_MAX),sum_hops(MY_INFINITY){}
    float min_cost;                 // for shortest path
    // predecessor;                 // for shortest path successor,
    uint64_t successor;             // for shortest path 
    float sum_distance;             // new
    uint64_t sum_hops;                           // new 
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

    bool at_KSPaths;             //vector<size_t> KSPaths_record;
    uint64_t predecessor;        // for shortest path, store the id of the graph tau
    uint64_t internel_id_of_g;   // internel id of orginal graph g

    float min_cost;              // for shortest path
    float sum_distance;          // for shortest path
    size_t sum_hops;             // for shortest path
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
 

/*bool is_loopless_path(graph_t& g, graph_tau& tau, size_t path_last_id_tau, size_t path_first_id_tau)//based on set
{
    size_t tmpId = path_last_id_tau;    
    vertex_iterator_tau vit_tau_tmp = tau.find_vertex(tmpId); 
    size_t path_len = vit_tau_tmp->property().sum_hops + 1;
    set<size_t> path_id_set_g;
    path_id_set_g.insert( vit_tau_tmp->property().internel_id_of_g );
    while (tmpId != path_first_id_tau)
    {
        tmpId = vit_tau_tmp->property().predecessor;  
        vit_tau_tmp = tau.find_vertex(tmpId);
        path_id_set_g.insert( vit_tau_tmp->property().internel_id_of_g );
    }
    return path_len == path_id_set_g.size();
}*/
bool is_loopless_path(graph_t& g, graph_tau& tau, size_t path_last_id_tau, size_t path_first_id_tau)//based on sorted vector
{
    size_t tmpId = path_last_id_tau;    
    vertex_iterator_tau vit_tau_tmp = tau.find_vertex(tmpId); 

    vector<size_t> path_id_set_g;  
    path_id_set_g.push_back( vit_tau_tmp->property().internel_id_of_g );
    while (tmpId != path_first_id_tau)
    {
        tmpId = vit_tau_tmp->property().predecessor;  
        vit_tau_tmp = tau.find_vertex(tmpId);
        path_id_set_g.push_back( vit_tau_tmp->property().internel_id_of_g );
    }
    sort(path_id_set_g.begin(), path_id_set_g.end());
    return adjacent_find(path_id_set_g.begin(), path_id_set_g.end()) == path_id_set_g.end();
}


void top_ksp_subFun(bool trueMinCost, size_t trueMinCost_Iter, size_t& curr_kValue, ofstream& myfile, graph_t& g, size_t src, size_t dest, size_t Kvalue, double max_phy_dist, size_t max_phy_hops, size_t min_phy_hops) // src and dest are exID
{

    uint64_t internel_src  = g.external_to_internel_id(to_string(src)); 
    uint64_t internel_dest = g.external_to_internel_id(to_string(dest)); 



    //// Now all processing is based on internel id first.
    //// (1) sssp sub-procedure  -->for vertex: update  v_vit->property().min_cost and v_vit->property().successor 
    priority_queue<pair_IntFlt, vector<pair_IntFlt>, min_comp_IntFlt> PQ;

    vertex_iterator dest_vit = g.find_vertex(internel_dest);  // note: here the source of sssp is internel_dest, 
    dest_vit->property().min_cost = 0;
    dest_vit->property().sum_hops = 0;     // new
    dest_vit->property().sum_distance = 0; // new
    PQ.push(pair_IntFlt(internel_dest,0));

    // vit->property().successor is used to construct sssp, where the ancestor of all nodes is internel_dest, 
    // by using vit->property().successor recursively, we can find the shortest path of any node to internel_dest
    while (!PQ.empty())   // sum_distance  sum_hops
    {
        size_t u = PQ.top().first; //id
        PQ.pop();

        vertex_iterator u_vit = g.find_vertex(u);
        for (edge_iterator eit = u_vit->in_edges_begin(); eit != u_vit->in_edges_end(); eit++)  // new for in edge
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);

            // for every  vertex, try relaxing the path
            if (trueMinCost)
            {
                float alt = u_vit->property().min_cost + eit->property().cost; // 
                if (alt < v_vit->property().min_cost) 
                {
                    v_vit->property().successor = u;  // new, ori is predecessor
                    v_vit->property().min_cost     = alt; 
                    v_vit->property().sum_hops     = u_vit->property().sum_hops + 1; 
                    v_vit->property().sum_distance = u_vit->property().sum_distance + eit->property().phy_dist; 
                    PQ.push(pair_IntFlt(v,alt));
                }
            }
            else
            {
                // min_cost  -- sum_hops    exchange
                int   alt = u_vit->property().sum_hops + 1; // 
                if (alt < v_vit->property().sum_hops) 
                {
                    v_vit->property().successor = u;  // new, ori is predecessor
                    v_vit->property().sum_hops     = alt; 
                    v_vit->property().min_cost     = u_vit->property().min_cost + eit->property().cost; 
                    v_vit->property().sum_distance = u_vit->property().sum_distance + eit->property().phy_dist; 
                    PQ.push(pair_IntFlt(v,alt));
                }
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

            float reduced_cost;
            if (trueMinCost)
            {
                reduced_cost= v_vit->property().min_cost - u_vit->property().min_cost + eit->property().cost;
            }
            else
            {
                reduced_cost = v_vit->property().sum_hops - u_vit->property().sum_hops + 1;  // min_cost  -- sum_hops    exchange
            }
            eit->property().reduced_cost = reduced_cost;
            u_vit->property().sorted_edges_of_vertex.push_back(pair_IntFlt(v,reduced_cost));     // pair: id (outedge), reduced_cost
        }
        sort(u_vit->property().sorted_edges_of_vertex.begin(), u_vit->property().sorted_edges_of_vertex.end(),
            [](pair_IntFlt vecInput_1, pair_IntFlt vecInput_2) {return (vecInput_1.second < vecInput_2.second);} );
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
    src_vit_tau->property().min_cost += g.find_vertex( src_vit_tau->property().internel_id_of_g )->property().min_cost;
    PQ_KSP_candidates_tau.push( pair_FltInt(src_vit_tau->property().min_cost, src_vit_tau->id()) );
    src_vit_tau->property().min_cost -= g.find_vertex( src_vit_tau->property().internel_id_of_g )->property().min_cost;


    size_t max_iter;
    if (trueMinCost)
    {
        max_iter = trueMinCost_Iter;
    }
    else
    {
        max_iter = Kvalue*100000;;
    }
    size_t k = curr_kValue;
    size_t iter = 0; 
    while (k<Kvalue && !PQ_KSP_candidates_tau.empty() && iter<max_iter)
    {
        // find the current shortest path p (may have loop)
        size_t candi_last_id = PQ_KSP_candidates_tau.top().second; 

        // X = X - {p}
        PQ_KSP_candidates_tau.pop();        

        // new
        uint64_t  cur_Gnode_atTau = tau.find_vertex(candi_last_id)->property().internel_id_of_g; 
        if (cur_Gnode_atTau != internel_dest)
        {
            vertex_iterator_tau dest_vit_tau = add_partialSP_totau(g, cur_Gnode_atTau, internel_dest, tau, candi_last_id, max_phy_dist, max_phy_hops);
            candi_last_id = dest_vit_tau->id();
        }
        
        // if p is loopless and within max_phy_dist and max_phy_hops, add it to the final KSP,  
        bool is_valid_candidate = (tau.find_vertex(candi_last_id)->property().internel_id_of_g == internel_dest); // new for the break in add_partial_candi 
        bool within_max_phy_dist = tau.find_vertex(candi_last_id)->property().sum_distance <= max_phy_dist;
        bool within_max_phy_hops = tau.find_vertex(candi_last_id)->property().sum_hops     <= max_phy_hops;
        bool largerthan_min_phy_hops = tau.find_vertex(candi_last_id)->property().sum_hops >= min_phy_hops;

        within_max_phy_dist = true; //Note!!!  for single-layer since no dist in this case.
        if ( is_valid_candidate && within_max_phy_dist && within_max_phy_hops && largerthan_min_phy_hops && is_loopless_path(g, tau, candi_last_id, internel_src_tau) )  
        {
            KSPaths_lastID_tau.push_back(candi_last_id);
            k++;
        }

        //note: even the current shortest path p may have loop, the new generated candidates based on p may have no loop.
        // find the deviation path pk_vkt_nodes, the top of pk_vkt_nodes is the deviation node
        stack<size_t> pk_vkt_nodes;
        size_t tmp_id = candi_last_id; // KSPaths_lastID_tau.back()
        vertex_iterator_tau dest_vit_tau = tau.find_vertex(tmp_id);
        dest_vit_tau->property().at_KSPaths = true;  
        while (tmp_id != internel_src_tau)  
        {
            tmp_id = dest_vit_tau->property().predecessor;

            pk_vkt_nodes.push(tmp_id);
            dest_vit_tau = tau.find_vertex(tmp_id);
            if (dest_vit_tau->property().at_KSPaths == true) 
            {
                break;
            }
            dest_vit_tau->property().at_KSPaths = true;  
        }

        // for each node in deviation path, try to find a candidate
        while (!pk_vkt_nodes.empty())  // for each node at pk_vkt
        {
            // find current deviation point: pk_top_id
            size_t pk_top_id = pk_vkt_nodes.top();  // tmp_id--> pk_top_id
            pk_vkt_nodes.pop();

            // get out if there is loop to save time
            if ( !is_loopless_path(g, tau, pk_top_id, internel_src_tau) )
            {
                break;
            }

            // find A(v), which is stored in neighborNodes_g (the first of each pair in the vector) and sorted based on reduced cost  
            vertex_iterator_tau pk_top_vit = tau.find_vertex(pk_top_id);   //vit_tau_tmp--> pk_top_vit
            size_t cur_deviation_id_g = pk_top_vit->property().internel_id_of_g;  // cur_deviation_id_g is the v in MPS
            vertex_iterator tmp_vit = g.find_vertex(cur_deviation_id_g);        
            vector<pair_IntFlt> neighborNodes_g = tmp_vit->property().sorted_edges_of_vertex;

            // (new version )find A_Tk_(v),  now further update it that  A_Tk_(v) contains out edges (sucsessors) and all predecessors up to internel_src_tau
            // this update is to reduce the number of loop paths
            vector<uint64_t> neighborNodes_tau; 
            // vector<uint64_t> neighborNodes_tau_at_g; //ori: 
            unordered_set<uint64_t> neighborNodes_tau_at_g; // to speed up the find operation later
            tau.find_vertex_out_neighborNodes(pk_top_id, neighborNodes_tau); // here only need  out neigbor, the
            for (size_t idx=0; idx<neighborNodes_tau.size(); ++idx)
            {
                vertex_iterator_tau vit_tau_tmptmp = tau.find_vertex(neighborNodes_tau[idx]);
                //neighborNodes_tau_at_g.push_back( vit_tau_tmptmp->property().internel_id_of_g );
                neighborNodes_tau_at_g.insert( vit_tau_tmptmp->property().internel_id_of_g );
            }
            tmp_id = pk_top_id;
            vertex_iterator_tau vit_tau_tmp = pk_top_vit;
            while (tmp_id != internel_src_tau)
            {
                tmp_id = vit_tau_tmp->property().predecessor;  
                vit_tau_tmp = tau.find_vertex(tmp_id);
                //neighborNodes_tau_at_g.push_back(vit_tau_tmp->property().internel_id_of_g);
                neighborNodes_tau_at_g.insert(vit_tau_tmp->property().internel_id_of_g);
            }


            // check the first arc in A(v) - A_Tk_(v) at Graph g:   neighborNodes_g - neighborNodes_tau_at_g             
            for (size_t idx=0; idx<neighborNodes_g.size(); ++idx)
            {
                //if ( !(find(neighborNodes_tau_at_g.begin(), neighborNodes_tau_at_g.end(), neighborNodes_g[idx].first) != neighborNodes_tau_at_g.end()) ) // try unordered_set
                if (  neighborNodes_tau_at_g.find(neighborNodes_g[idx].first) == neighborNodes_tau_at_g.end()  ) // --> similar or even slow in small cases.
                {
                    // v and x is the cur_deviation_id_g and neighborNodes_g[idx].first, respectively, in MPS alg.
                    // now add point x and edge v->x in tau  vertex_iterator    
                    vertex_iterator_tau dest_vit_tau = tau.add_vertex(); // this is x in tau

                    edge_iterator eit_g;
                    bool find_result = g.find_out_edge_2id(cur_deviation_id_g, neighborNodes_g[idx].first, eit_g); // for eit_g->property().cost and eit_g->property().phy_dist
                    assert(find_result);

                    dest_vit_tau->property().at_KSPaths = false;//dest_vit_tau->property().KSPaths_record.push_back(1); // this node is at the 1st shortest path.
                    dest_vit_tau->property().predecessor = pk_top_id; // pk_top_id = pk_vkt_nodes.top();
                    dest_vit_tau->property().internel_id_of_g = neighborNodes_g[idx].first;   
                    dest_vit_tau->property().min_cost     = pk_top_vit->property().min_cost + eit_g->property().cost; //   
                    dest_vit_tau->property().sum_distance = pk_top_vit->property().sum_distance + eit_g->property().phy_dist; //   
                    dest_vit_tau->property().sum_hops     = pk_top_vit->property().sum_hops + 1;

                    uint64_t tau_tmp_id = dest_vit_tau->id();  //  
                    edge_iterator_tau eit_tau;
                    tau.add_edge(pk_top_id,tau_tmp_id,eit_tau); // note: put all info at vertex, see dest_vit_tau  (v,x)

   
                    dest_vit_tau->property().min_cost += g.find_vertex( dest_vit_tau->property().internel_id_of_g )->property().min_cost;
                    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().min_cost, dest_vit_tau->id()) );
                    dest_vit_tau->property().min_cost -= g.find_vertex( dest_vit_tau->property().internel_id_of_g )->property().min_cost;
                    //dest_vit_tau->property().sum_distance -= vit_g->property().sum_distance;
                    //dest_vit_tau->property().sum_hops     -= vit_g->property().sum_hops;

                    break;
                } 
            } 
        }//while (!pk_vkt_nodes.empty())  // for each node at pk_vkt
        
        iter++;

    } //while (k<Kvalue && !PQ_KSP_candidates_tau.empty() && iter<max_iter)
    curr_kValue = k;




    //#ifdef TOPKSP_PRINTOUT
    // output_top_ksp: 
    //cout<<"The top "<<Kvalue<<" shortest loopless paths";
    cout<<"Output: trueMinCost is "<< trueMinCost <<endl;;
    cout<<"--the running iteration (do PQ pop) is "<<iter<<endl;
    cout<<"--the found final path number is "<<KSPaths_lastID_tau.size()<<endl;

    if (KSPaths_lastID_tau.size()==0)
    {
        cout<<"cannot find any path, maybe we have too much constraints!"<<endl;
        assert(false);
    }

    
    for (size_t idx = 0; idx < KSPaths_lastID_tau.size(); ++idx)
    {
        size_t tmp_k = idx + 1;
        uint64_t curr_id_tau = KSPaths_lastID_tau[idx];

        vertex_iterator_tau vit_tau = tau.find_vertex(curr_id_tau);
        /*
        cout<<endl;
        cout<<"The "<<tmp_k<<" result with metrics:   ";
        cout<<"min_cost is "    <<vit_tau->property().min_cost<<"; ";
        cout<<"sum_distance is "<<vit_tau->property().sum_distance<<"; ";
        cout<<"sum_hops is "    <<vit_tau->property().sum_hops<<"."<<endl;
        */

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

        //cout<<"The "<<tmp_k<<" result with path nodes: ";
        
        myfile<<" AggSrc,AggDst,Path"<<endl;
        myfile<<"0,2,";
        while (!curr_path.empty())
        {
            int tmp_node = curr_path.top();
            myfile<< tmp_node << "|";
            //cout<<"-->"<<tmp_node;
            curr_path.pop();
        }
        myfile<<",\n";

        if ( !is_loopless_path(g, tau, KSPaths_lastID_tau[idx], internel_src_tau) )
        {
            cout<<endl<<"wrong, this path has loop"<<endl;
            assert(false);
        }
        //cout<<endl;
    }
    

    /*
    if (Kvalue!=KSPaths_lastID_tau.size())
    {
        cout<<"Warning: the input k= "<<Kvalue<<" is too large, cannot find the required loopless paths with current constraints!"<<endl;
        cout<<"the largest possible k is "<<KSPaths_lastID_tau.size()<<endl;
    }
    */
    //#endif

    //cout<<"(the acutal number of running iteration is "<<iter<<"):"<<endl;

}

void top_ksp(graph_t& g, size_t src, size_t dest, size_t Kvalue, double max_phy_dist, size_t max_phy_hops, size_t min_phy_hops, gBenchPerf_event & perf, int perf_group)
{

    perf.open(perf_group);
    perf.start(perf_group);
#ifdef SIM
    SIM_BEGIN(true);  
#endif


    ofstream myfile ("outputSingleLayer/combineMethods_single_0_2.csv"); // update the above number
    size_t trueMinCost_Iter = 5*1e6; // 1e6   (0-334:  1e7->152s)

    bool trueMinCost        = true;
    size_t curr_kValue      = 0;


    // first operate trueMinCost with trueMinCost_Iter iterations
    top_ksp_subFun(trueMinCost, trueMinCost_Iter, curr_kValue, myfile, g, src, dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops);

    // second, find the remaining (Kvalue - curr_kValue) paths by using reduced hops
    reset_graph(g);
    trueMinCost = false;
    top_ksp_subFun(trueMinCost, trueMinCost_Iter, curr_kValue, myfile, g, src, dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops);

    myfile.close();

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
    arg.add_arg("min_phy_hops","1","root/min_phy_hops");
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

    size_t src, dest, Kvalue, max_phy_hops,min_phy_hops, threadnum;
    double max_phy_dist;
    arg.get_value("src",src);
    arg.get_value("dest",dest);
    arg.get_value("Kvalue",Kvalue);
    arg.get_value("max_phy_dist",max_phy_dist);
    arg.get_value("max_phy_hops",max_phy_hops);
    arg.get_value("min_phy_hops",min_phy_hops);
    arg.get_value("threadnum",threadnum);



    graph_t graph;
    cout<<"loading data... \n";
    t1 = timer::get_usec();
    int AGG_pair_num, total_vertex;





    srand (1223);
    //Kvalue = 5;  // 5   50    500  
    AGG_pair_num = 1;  // 5    50    100 

    max_phy_dist = FLT_MAX; // Note, since this is no distance for Single-layer.
    max_phy_hops = 9; //9
    min_phy_hops = 4; //4
    
    //  edge  vertex
    total_vertex = 400;  // 400   4000   10000 
    string vfile = path + "/vertex400.csv";  // vertex400   vertex4000    vertex10000
    string efile = path + "/vertex400dim8_edge.csv";       
    //  vertex400dim8_edge     vertex400dim12_edge       vertex400dim16_edge
    //  vertex4000dim8_edge    vertex4000dim12_edge      vertex4000dim16_edge
    //  vertex10000dim8_edge   vertex10000dim12_edge     vertex10000dim16_edge






// #ifndef EDGES_ONLY    
//     if (graph.load_csv_vertices(vfile, true, separator, 0) == -1)
//         return -1;
//     if (graph.load_csv_edges(efile, true, separator, 0, 1) == -1) 
//         return -1;
// #else
//     if (graph.load_csv_edges(efile, true, separator, 0, 1) == -1)
//         return -1;
// #endif









    vfile = path + "/Node.csv";  
    efile = path + "/Link.csv";    
    efile = path + "/LinkUndirected.csv";  // add same link two times for undirected graph
    if (graph.load_csv_vertices(vfile, true, ",", 0) == -1)
        return -1;
    //if (graph.load_csv_edges(efile, true, "," , 0, 1) == -1) 
    if (graph.load_csv_edges(efile, true, ",", 1, 2,false, NULL,3,-1) == -1) 
        return -1;






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
        {

            cout<<"max_phy_dist is " << max_phy_dist <<endl;
            cout<<"max_phy_hops is " << max_phy_hops <<endl;
            cout<<"min_phy_hops is " << min_phy_hops <<endl;

            for (int idx=0; idx<AGG_pair_num; idx++) //
            {
               // int min, max;
                //min = 0;
                //max = 19;
                //src = rand()%(max-min + 1) + min;
                
               // min  = total_vertex-20;
               // max  = total_vertex-1;                
               // dest = rand()%(max-min + 1) + min;


                //src  = 0;
                //dest = 1;



                #ifdef TOPKSP_PRINTOUT
                cout<<endl<<endl<<"Now output the top "<<Kvalue<<" Shortest paths bwteen nodes: "<<src<<"-->"<<dest;
                cout<<"(run "<<idx+1<<" out of "<<AGG_pair_num<<")"<<endl;
                #endif

                //src = 46;
                //dest = 252;
                top_ksp(graph, src , dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops, perf, i);  
                //reset_graph(graph);
            }

        //top_ksp(graph, src , dest, Kvalue, max_phy_dist, max_phy_hops, perf, i);  

        }
        else
        {
            //parallel_sssp(graph, root, threadnum, perf_multi, i);
        }


        t2 = timer::get_usec();
        elapse_time += t2-t1;
        if ((i+1)<run_num) reset_graph(graph);
    }
#ifndef ENABLE_VERIFY
    cout<<"AGG_pair_num is "<< AGG_pair_num <<"; Kvalue is "<<Kvalue<<endl;
    cout<<"== Total ruuning time: "<<elapse_time/run_num<<" sec\n";
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

