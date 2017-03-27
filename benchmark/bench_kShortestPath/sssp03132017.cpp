//====== Graph Benchmark Suites ======//
//========== Shortest Path ===========//
// 
// Single-source shortest path
// 
// Usage: ./sssp    --dataset <dataset path> 
//                  --root <root vertex id> 
//                  --target <targegt vertex id>

#include "common.h"
#include "def.h"
#include "openG.h"
#include "omp.h"
#include <fstream>
#include <queue>
#include <cfloat>     
#include <cassert>
#include <stack>
#include <cstdlib>
#include <unordered_set>

using namespace std;

#ifdef HMC
#include "HMC.h"
#endif

#ifdef SIM
#include "SIM.h"
#endif

//#define VERIFY_RESULTS     // if compare serial results with parallel results
#define MY_INFINITY  0xfff0  // for unsigned i.e., 65520

size_t beginiter = 0;
size_t enditer = 0;

typedef pair<size_t,float> pair_IntFlt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<float,size_t> pair_FltInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<size_t,size_t> pair_IntInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //

class vertex_property // new
{
public:  // sum_distance(0),sum_hops(0){}
    vertex_property():min_cost(FLT_MAX),successor(MY_INFINITY),sum_distance(FLT_MAX),sum_hops(MY_INFINITY),weight(FLT_MAX),occurrence(0){}
    float min_cost;                 // for shortest path
    // predecessor;                 // for shortest path successor,
    uint64_t successor;             // for shortest path 
    float sum_distance;             // new
    uint64_t sum_hops;              // new 
	float weight;                   // new
	uint64_t occurrence;            // new
	
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
		vit->property().weight       = FLT_MAX;
		vit->property().occurrence   = 0;
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
    vertex_property_tau():at_KSPaths(false),predecessor(MY_INFINITY),internel_id_of_g(MY_INFINITY),min_cost(0),sum_distance(0),sum_hops(0),weight(0),occurrence(0){}

    bool at_KSPaths;             //vector<size_t> KSPaths_record;
    uint64_t predecessor;        // for shortest path, store the id of the graph tau
    uint64_t internel_id_of_g;   // internel id of orginal graph g

    float min_cost;              // for shortest path
    float sum_distance;          // for shortest path
    size_t sum_hops;             // for shortest path
	float weight;                // for shortest path
	uint64_t occurrence;         // counting the occurrence of a node
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
#endif


vertex_iterator_tau add_partialSP_totau(graph_t& g, size_t src_id_g, size_t dest_id_g, graph_tau& tau, size_t start_id_tau,   \
double max_phy_dist, size_t max_phy_hops, size_t alpha)  
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
		dest_vit_tau->property().occurrence   = dest_vit_g->property().occurrence;
		
		unsigned int alpha_occur = (dest_vit_tau->property().occurrence > 7000) ? int(alpha*10) : int(alpha*10);
		dest_vit_tau->property().weight       = dest_vit_tau->property().min_cost + dest_vit_tau->property().sum_hops*alpha + dest_vit_tau->property().occurrence*alpha_occur;  //weight
		

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


void top_ksp_subFun(bool trueMinCost, size_t trueMinCost_Iter, size_t& curr_kValue, ofstream& myfile,   \
graph_t& g,  size_t src, size_t dest, size_t Kvalue, double max_phy_dist, size_t max_phy_hops,          \
size_t min_phy_hops, size_t alpha, unordered_set<string>& paths) // src and dest are exID
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
	dest_vit->property().weight = 0;       // new
	//dest_vit->property().occurrence = dest_vit->property().occurrence;       // new
	
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
                unsigned int alt = u_vit->property().sum_hops + 1; // 
				float alt_cost = u_vit->property().min_cost + eit->property().cost; //  
				unsigned int occur = v_vit->property().occurrence;
				
				unsigned int alpha_occur = (occur > 7000) ? int(alpha*10) : int(alpha*10);
				unsigned int weight = alt*alpha + alt_cost + occur*alpha_occur;
				
                //if (alt < v_vit->property().sum_hops) 
				if (weight < v_vit->property().weight) 
                {
                    v_vit->property().successor = u;  // new, ori is predecessor
                    v_vit->property().sum_hops     = alt; 
                    v_vit->property().min_cost     = u_vit->property().min_cost + eit->property().cost; 
                    v_vit->property().sum_distance = u_vit->property().sum_distance + eit->property().phy_dist; 
					v_vit->property().weight = weight;
					
                    PQ.push(pair_IntFlt(v,weight));
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
                //reduced_cost = v_vit->property().sum_hops - u_vit->property().sum_hops + 1;  // min_cost  -- sum_hops    exchange
				unsigned int hops_gap  = v_vit->property().sum_hops - u_vit->property().sum_hops + 1; 
				//generalized reduced_cost = edge_weight + edge_cost/alpha + alpha*(num_hops=1) + dest_v.occurrence *(alpha*10.0)
				
				unsigned int alpha_occur = (v_vit->property().occurrence > 7000) ? int(alpha*10.0) : int(alpha*10.0);
				reduced_cost = v_vit->property().min_cost - u_vit->property().min_cost + eit->property().cost + alpha*(hops_gap) + v_vit->property().occurrence*alpha_occur; 
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
	src_vit_tau->property().weight = 0;
	src_vit_tau->property().occurrence = g.find_vertex(internel_src)->property().occurrence;
	
    // construct the first shortest path constructed in tau
    uint64_t internel_src_tau = src_vit_tau->id();
    src_vit_tau->property().min_cost += g.find_vertex( src_vit_tau->property().internel_id_of_g )->property().min_cost;
	src_vit_tau->property().weight += g.find_vertex( src_vit_tau->property().internel_id_of_g )->property().weight;
    PQ_KSP_candidates_tau.push( pair_FltInt(src_vit_tau->property().weight, src_vit_tau->id()) );   //????
    src_vit_tau->property().min_cost -= g.find_vertex( src_vit_tau->property().internel_id_of_g )->property().min_cost;
	src_vit_tau->property().weight -= g.find_vertex( src_vit_tau->property().internel_id_of_g )->property().weight;


    size_t max_iter;
    if (trueMinCost)
    {
        max_iter = trueMinCost_Iter;
    }
    else
    {
        //max_iter = Kvalue*100000;;
		max_iter = trueMinCost_Iter;
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
            vertex_iterator_tau dest_vit_tau = add_partialSP_totau(g, cur_Gnode_atTau, internel_dest, tau, candi_last_id, max_phy_dist, max_phy_hops, alpha);
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
			//Obtain the key for the path
			uint64_t candi_id_tau = candi_last_id;  
			bool overused = 0;
			string key = "";
			vertex_iterator_tau vit_tau = tau.find_vertex(candi_id_tau);
			vector<vertex_iterator> iters_candi_path;

			uint64_t candi_id_g = vit_tau->property().internel_id_of_g; 
			uint64_t candi_exID = atoi(g.internal_to_externel_id(candi_id_g).c_str());
			key += to_string(candi_exID);
			do
			{
				candi_id_tau = vit_tau->property().predecessor;
				vit_tau = tau.find_vertex(candi_id_tau);
				candi_id_g = vit_tau->property().internel_id_of_g; 
				candi_exID = atoi(g.internal_to_externel_id(candi_id_g).c_str());
				
				//Obtain the key for the path
				key += to_string(candi_exID);
				//Obtain iterators for all vertexs in the path
				vertex_iterator v_iter_candi_path = g.find_vertex(candi_exID);
				iters_candi_path.push_back(v_iter_candi_path);
				//check if there exits an overused node
				if (v_iter_candi_path->property().occurrence > 40000)
				{
					//cout<<v_iter_candi_path->property().occurrence<<"?\n";
					overused = 1;
				}
				
			} while (candi_id_tau != internel_src_tau);
			
			if (paths.find(key) == paths.end() and !overused)
			{
				KSPaths_lastID_tau.push_back(candi_last_id);
				k++;
				//insert a new key
				paths.insert(key);
				vector<vertex_iterator>::iterator v_iter;
				//update occurrence of each vertex
				for (v_iter = iters_candi_path.begin(); v_iter != iters_candi_path.end(); v_iter++)
				{
					(*v_iter)->property().occurrence += 1;
				}
			}
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
                    
                    unsigned int alpha_occur = (dest_vit_tau->property().occurrence > 7000) ? int(alpha*10.0) : int(alpha*10.0);
                    dest_vit_tau->property().weight     = dest_vit_tau->property().min_cost + dest_vit_tau->property().sum_hops*alpha + dest_vit_tau->property().occurrence*alpha_occur ;

                    uint64_t tau_tmp_id = dest_vit_tau->id();  //  
                    edge_iterator_tau eit_tau;
                    tau.add_edge(pk_top_id,tau_tmp_id,eit_tau); // note: put all info at vertex, see dest_vit_tau  (v,x)

   
                    dest_vit_tau->property().min_cost += g.find_vertex( dest_vit_tau->property().internel_id_of_g )->property().min_cost;
                    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().weight, dest_vit_tau->id()) );    //????
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

    cout<<"Output: trueMinCost is "<< trueMinCost <<endl;;
    cout<<"--the running iteration (do PQ pop) is "<<iter<<endl;
    cout<<"--the found final path number is "<<KSPaths_lastID_tau.size()<<endl;

    /*if (KSPaths_lastID_tau.size()==0)
    {
        cout<<"cannot find any path, maybe we have too much constraints!"<<endl;
        assert(false);
    
    }*/

    for (size_t idx = 0; idx < KSPaths_lastID_tau.size(); ++idx)
    {
        //size_t tmp_k = idx + 1;
        uint64_t curr_id_tau = KSPaths_lastID_tau[idx];

        vertex_iterator_tau vit_tau = tau.find_vertex(curr_id_tau);

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
        //myfile<<"0,2,";
        myfile<<src<<","<<dest<<",";
        
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
 
 
void top_ksp(size_t trueMinCost_Iter, ofstream& myfile, graph_t& g, size_t src, size_t dest,  \
size_t Kvalue, double max_phy_dist, size_t max_phy_hops, size_t min_phy_hops, 
gBenchPerf_event & perf, int perf_group )
{

    perf.open(perf_group);
    perf.start(perf_group);
#ifdef SIM
    SIM_BEGIN(true);  
#endif

    double sum_cost     = 0;
    size_t number_edges = 0;
    vertex_iterator vit;
    for (vit=g.vertices_begin(); vit!=g.vertices_end(); vit++)
    {
        for (edge_iterator eit = vit->in_edges_begin(); eit != vit->in_edges_end(); eit++)  // new for in edge
        {
            if ( abs(FLT_MAX - eit->property().cost) > FLT_MAX/2 )
            {
                sum_cost += eit->property().cost;
                ++number_edges;
            }
        }
    }
    unsigned int rate = int(sum_cost/number_edges);
    
    bool trueMinCost        = true;
    size_t curr_kValue      = 0;
    unordered_set<string> paths;

    // first operate trueMinCost with trueMinCost_Iter iterations
    //top_ksp_subFun(trueMinCost, trueMinCost_Iter, curr_kValue, myfile, g, src, dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops, 1, paths);

    // second, find the remaining (Kvalue - curr_kValue) paths by using reduced hops
    int alpha = 1;
    while (curr_kValue < Kvalue)
    {
        reset_graph(g);
        trueMinCost = false;
        alpha *= 8; //int(math.sqrt(rate));
        top_ksp_subFun(trueMinCost, trueMinCost_Iter, curr_kValue, myfile, g, src, dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops, alpha, paths);
    }
    //myfile.close();
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
    arg.add_arg("min_phy_hops","1","root/min_phy_hops");
}
//==============================================================//


//Load AGG_pairs
vector<vector<int> > loadAggPairs(void)
{
    vector<vector<int> > res;

   // new data
    int raw[21][3] = {{0,  1,  18716},
                      {0,  2,  30006},
                      {0 , 3 , 16964},
                      {0     , 334   , 42659 },
                      {0     , 394   , 42834 },
                      {0     , 557   , 42983 },

                      {1     , 2     , 31163 },
                      {1     , 3     , 25949 },
                      {1     , 334   , 46181 },
                      {1     , 394   , 46280 },

                      {2     , 557   , 42043 },
                      {1     , 557   , 46006 },
                      {2     , 3     , 21733 },
                      {2     , 334   , 41844 },
                      {2     , 394   , 41884 },

                      {3     , 334   , 43992 },
                      {3     , 394   , 44119 },
                      {3     , 557   , 44324 },
                      {334   , 394   , 29481 },
                      {334   , 557   , 29358 },
                      {394   , 557   , 28410 }};
    
    for (int i=0; i<21; i++)
    {
        vector<int> pair;
        for (int j=0; j<3; j++)
            pair.push_back(raw[i][j]);
        res.push_back(pair);
    }
    return res;

}


//Initialize a graph
void graphInit(graph_t &graph, string vfile, string efile, size_t src, size_t dest, string separator)
{
    cout<<"loading data... \n";
    double t1 = timer::get_usec();
    
    if (graph.load_csv_vertices(vfile, true, ",", 0) == -1)
        return;
    if (graph.load_csv_edges(efile, true, ",", 1, 2,false, NULL,3,-1) == -1) 
        return;

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


double serialTest(size_t true_min_iter, string vfile, string efile, size_t src, size_t dest, string separator  \
, int max_phy_dist, int max_phy_hops, size_t min_phy_hops,  vector<vector<int> >tests, gBenchPerf_event perf)
{
    graph_t *graph;  
    double t1 = timer::get_usec();
    int p_src, p_dest, k_val;

    ofstream res_fstream;
    res_fstream.open("./outputSingleLayer/combineMethods.csv", ofstream::trunc);
     
    for (unsigned int i=0; i<tests.size(); i++)
    {
        graph = new graph_t();
        graphInit(*graph, vfile, efile, src, dest, separator);
        p_src = tests[i][0]; 
        p_dest = tests[i][1];
        k_val = tests[i][2];
        
        //res_fstream.open("./test_res/parallel_res_"+to_string(idx), ofstream::trunc);
        top_ksp(true_min_iter, res_fstream, *graph, p_src , p_dest, k_val, max_phy_dist, max_phy_hops, min_phy_hops, perf, i);  
        delete graph;
    }

    res_fstream.close();
    double t2 = timer::get_usec();
    return t2-t1;
}



double parallelTest(size_t true_min_iter, string vfile, string efile, size_t src, size_t dest, string separator,    \
int max_phy_dist, int max_phy_hops, size_t min_phy_hops, vector<vector<int> >tests, gBenchPerf_event perf)
{
    graph_t *graph;     
    double t1 = timer::get_usec();

    omp_set_num_threads(tests.size());
    #pragma omp parallel private(graph)
    {
        int idx = omp_get_thread_num();
        int p_src = tests[idx][0];
        int p_dest = tests[idx][1];
        int k_val = tests[idx][2];
        
        ofstream res_fstream;
        res_fstream.open("outputSingleLayer/combineMethods_parallel"+to_string(idx)+".csv");

        graph = new graph_t();
        graphInit(*graph, vfile, efile, src, dest, separator);

        top_ksp(true_min_iter, res_fstream, *graph, p_src , p_dest, k_val, max_phy_dist, max_phy_hops, min_phy_hops, perf, 0);  
        delete graph;

        res_fstream.close();
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
    size_t src, dest, Kvalue, max_phy_hops, min_phy_hops;
    size_t threadnum, AGG_pair_num;
    double max_phy_dist;
    string vfile, efile, aggfile;

    arg.get_value("dataset",path);
    arg.get_value("separator",separator);
    arg.get_value("src",src);
    arg.get_value("dest",dest);
    arg.get_value("Kvalue",Kvalue);
    arg.get_value("max_phy_dist",max_phy_dist);
    arg.get_value("max_phy_hops",max_phy_hops);
    arg.get_value("min_phy_hops",min_phy_hops);
    arg.get_value("threadnum",threadnum);

    
    AGG_pair_num = 2;   // 5    50    100 
    Kvalue = 30;         // 5   50    500  
    max_phy_dist = FLT_MAX;
    max_phy_hops = 9;
    min_phy_hops = 4;
    size_t true_min_iter =  1e5; //1e5; //5*1e6; // 1e6   (0-334:  1e7->152s)
     
    //int total_vertex = 400;  // 400   4000   10000 
    vfile = path + "/Node.csv";  // vertex400   vertex4000    vertex10000
    efile = path + "/LinkUndirected.csv";       
    aggfile = path + "/AggPair.csv";
    
    
    cout<<"...\n";
    gBenchPerf_multi perf_multi(threadnum, perf);

    //Generate tests
    //vector<vector<int> > tests = genTests(total_vertex, AGG_pair_num);
    vector<vector<int> > tests = loadAggPairs();
    AGG_pair_num = tests.size();
    
    double s_time, p_time;

    cout<<"AGG_pair_num is "<< AGG_pair_num <<"; Kvalue is "<<Kvalue<<endl;
    cout<<"max_phy_dist is " << max_phy_dist <<endl;
    cout<<"max_phy_hops is " << max_phy_hops <<endl;
    cout<<"min_phy_hops is " << min_phy_hops <<endl<<endl;
    
    //Serial run
    cout<<"Start running serial test"<<endl;
    s_time = serialTest(true_min_iter, vfile, efile, src, dest, separator, max_phy_dist, max_phy_hops, min_phy_hops, tests, perf);

    #ifndef ENABLE_VERIFY
        cout<<"== Total ruuning time: "<<s_time<<" sec\n";
        cout<<"==================================================================\n"<<endl;
    #endif    

    //Parallel run
    cout<<"Start running openMP test"<<endl;
    p_time = parallelTest(true_min_iter, vfile, efile, src, dest, separator, max_phy_dist, max_phy_hops, min_phy_hops, tests, perf);
    
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

