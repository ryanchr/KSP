//====== Graph Benchmark Suites ======//
//========== Shortest Path ===========//
// 
// Single source-dest pair shortest path with multi-constraints
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
#define TOPKSP_DEBUG

size_t beginiter = 0;
size_t enditer = 0;

typedef pair<size_t,float> pair_IntFlt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<float,size_t> pair_FltInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //
typedef pair<size_t,size_t> pair_IntInt; //  //typedef pair<size_t,size_t> pair_IntFlt; //


class vertex_property // new
{
public:  // sum_distance(0),sum_hops(0){}
    vertex_property():min_cost(FLT_MAX),successor(MY_INFINITY),sum_distance(FLT_MAX),sum_hops(MY_INFINITY)  \
    ,weight(FLT_MAX),occurrence(0),visited(0){}
    float min_cost;                 // for shortest path
    // predecessor;                 // for shortest path successor,
    uint64_t successor;             // for shortest path 
    float sum_distance;             // new
    uint64_t sum_hops;              // new 
	float weight;                   // new
	uint64_t occurrence;            // new
    bool visited;                   // new
	
    vector<pair_IntFlt> sorted_edges_of_vertex;  // pair: id (outedge), reduced_cost
};

class edge_property // new
{
public:
    edge_property():cost(FLT_MAX),phy_dist(FLT_MAX),max_bw(FLT_MAX),up_bw_residue(FLT_MAX),
                    down_bw_residue(FLT_MAX),reduced_cost(FLT_MAX){} // new:  Note: 

    float cost; // note: here cost means cost
    float phy_dist; //new
    float max_bw;
    float up_bw_residue;
    float down_bw_residue;
    float reduced_cost; 
};

typedef openG::extGraph<vertex_property, edge_property> graph_t;
typedef graph_t::vertex_iterator    vertex_iterator;
typedef graph_t::edge_iterator      edge_iterator;


void reset_graph(graph_t & g, bool reset_occurrence)
{
    vertex_iterator vit;
    for (vit=g.vertices_begin(); vit!=g.vertices_end(); vit++)
    {
        vit->property().min_cost     = FLT_MAX;
        vit->property().successor    = MY_INFINITY;
        vit->property().sum_distance = FLT_MAX;
        vit->property().sum_hops     = MY_INFINITY;
		vit->property().weight       = FLT_MAX;
        vit->property().visited      = 0;
        if (reset_occurrence)
		  vit->property().occurrence   = 0;
        vit->property().sorted_edges_of_vertex.clear();

        for (edge_iterator eit = vit->in_edges_begin(); eit != vit->in_edges_end(); eit++)  // new for in edge
        {
            eit->property().reduced_cost = FLT_MAX;
        }
    }
}

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

    void output_shorest_path(graph_t& g, size_t src, size_t dest) //for test, src and dest is exID, dest is the input dest of spmc func.
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
            //cout<<"-->"<<curr_exID<<endl;
        } while (curr_id_g != internel_dest);
        cout<<endl;
    }

    void printVec(vector<size_t> in)
    {
        for (auto iter = in.begin(); iter != in.end(); iter++)
            cout<<*iter<<",";
        cout<<endl;
    }
#endif


//calculate the shortest path of any node to internel_dest
void sssp(graph_t& g,  size_t src, size_t dest, size_t alpha)
{
    //uint64_t internel_src  = g.external_to_internel_id(to_string(src)); 
    uint64_t internel_dest = g.external_to_internel_id(to_string(dest)); 
    priority_queue<pair_IntFlt, vector<pair_IntFlt>, min_comp_IntFlt> PQ;

    //// Now all processing is based on internel id first.
    //// (1) sssp sub-procedure  -->for vertex: update  v_vit->property().min_cost and v_vit->property().successor 
    vertex_iterator dest_vit = g.find_vertex(internel_dest);  // note: here the source of sssp is internel_dest, 
    dest_vit->property().min_cost = 0;
    dest_vit->property().sum_hops = 0;     // new
    dest_vit->property().sum_distance = 0; // new
	dest_vit->property().weight = 0;       // new
    dest_vit->property().visited = 0;
	//dest_vit->property().occurrence = dest_vit->property().occurrence;       // new
	
    PQ.push(pair_IntFlt(internel_dest,0));

    // vit->property().successor is used to construct sssp, where the ancestor of all nodes is internel_dest, 
    // by using vit->property().successor recursively, we can find the shortest path of any node to internel_dest
    while (!PQ.empty())   // sum_distance  sum_hops
    {
        size_t u = PQ.top().first; //id
        PQ.pop();

        vertex_iterator u_vit = g.find_vertex(u);
        //in_edges: u is the source, v is the destination
        for (edge_iterator eit = u_vit->in_edges_begin(); eit != u_vit->in_edges_end(); eit++)  // new for in edge
        {
            size_t v = eit->target();

            //if (u == 15 && src == 0 && dest == 18)
            //    cout<<"Out edges of node 15: "<<v<<",";

            vertex_iterator v_vit = g.find_vertex(v);

            // for every  vertex, try relaxing the path
            // min_cost  -- sum_hops    exchange
            unsigned int alt = u_vit->property().sum_hops + 1; // 
			//float alt_cost = u_vit->property().min_cost + eit->property().cost; //  
            float alt_weight = u_vit->property().weight + eit->property().cost;

			unsigned int occur = v_vit->property().occurrence;
			
			unsigned int alpha_occur = (occur > 7000) ? int(alpha*10) : int(alpha*10);
			unsigned int weight = alt_weight + occur*alpha_occur;
			
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
void updateReducedCost(graph_t& g, size_t alpha)
{
    for (vertex_iterator u_vit=g.vertices_begin(); u_vit!=g.vertices_end(); u_vit++) // for each vertex u
    {
        for (edge_iterator eit = u_vit->edges_begin(); eit != u_vit->edges_end(); eit++) // for each outedge u->v from vertex u
        {
            size_t v = eit->target();
            vertex_iterator v_vit = g.find_vertex(v);
           
            //reduced_cost = v_vit->property().sum_hops - u_vit->property().sum_hops + 1;  // min_cost  -- sum_hops    exchange
			//unsigned int hops_gap  = v_vit->property().sum_hops - u_vit->property().sum_hops + 1; 
			//generalized reduced_cost = edge_weight + edge_cost/alpha + alpha*(num_hops=1) + dest_v.occurrence *(alpha*10.0)
			//int alpha_occur = (v_vit->property().occurrence > 7000) ? int(alpha*10.0) : int(alpha*10.0);
			
            //float reduced_cost = v_vit->property().min_cost - u_vit->property().min_cost + eit->property().cost + v_vit->property().occurrence*alpha_occur; 
            float reduced_weight = v_vit->property().weight - u_vit->property().weight + eit->property().cost; 
            
            eit->property().reduced_cost = reduced_weight;
            u_vit->property().sorted_edges_of_vertex.push_back(pair_IntFlt(v,reduced_weight));     // pair: id (outedge), reduced_cost
        }

        sort(u_vit->property().sorted_edges_of_vertex.begin(), u_vit->property().sorted_edges_of_vertex.end(),
            [](pair_IntFlt vecInput_1, pair_IntFlt vecInput_2) {return (vecInput_1.second < vecInput_2.second);} );
    }
}


//Return the shortest path for a single src-dest pair
void getShortestPath(graph_t& g, size_t src, size_t dest, size_t min_bw, vector<size_t> &shorest_path)
{
    uint64_t internel_src = g.external_to_internel_id(to_string(src));
    uint64_t internel_dest = g.external_to_internel_id(to_string(dest));
    
    vertex_iterator v_cur_node = g.find_vertex(internel_src);
    size_t cur_node = internel_src;

    int num_iter = 0;
    cout<<"Src: "<<src<<", Dest: "<<dest<<endl;
    while (cur_node != dest and num_iter < 25)
    {
        shorest_path.push_back(cur_node);
        vector<pair_IntFlt> &neighbor_edges = v_cur_node->property().sorted_edges_of_vertex;
        auto edge_iter = neighbor_edges.begin();

        for (; edge_iter != neighbor_edges.end(); ++edge_iter)
        {
            edge_iterator cur_edge;
            size_t cur_v_target = edge_iter->first;
            g.find_out_edge_2id(cur_node, cur_v_target, cur_edge);

            if (cur_edge->property().down_bw_residue >= min_bw)
            {
                cur_node = cur_v_target;
                break;
            }
        }

        if(edge_iter == neighbor_edges.end())
            cout<<"Cannot find a shortest path";
            assert(0);

        v_cur_node = g.find_vertex(cur_node);
        num_iter++;
    }

    shorest_path.push_back(internel_dest);
}


//back tracking search
void backTrackingSearch(graph_t& g, unordered_set<int>& invalid_nodes, size_t src, size_t dest
                       ,size_t min_bw, vector<size_t>& stack_nodes, bool &findPath)
{
    //Check if tmp_node to dest is another path sharing no edge with shortest path
    bool debug = 0; //(src == 15 && dest == 8);

    if (debug)
        printVec(stack_nodes);

    size_t tmp_node = stack_nodes.back();
    bool no_share_nodes = 1;
    bool has_bw = 1;
    vector<int> nodes_to_add;

    //int num_iter = 0;

    while (tmp_node != dest ) //&& num_iter < 30)
    {
        //num_iter++;

        vertex_iterator v_tmp_node = g.find_vertex(tmp_node);

        if (invalid_nodes.find(int(tmp_node)) != invalid_nodes.end())
        {
            v_tmp_node->property().occurrence += 1;
            no_share_nodes = 0;
            break;
        }

        size_t pre_node = tmp_node;
        edge_iterator edge_iter;

        tmp_node = v_tmp_node->property().successor;
        nodes_to_add.push_back(tmp_node);

        bool find_edge = g.find_out_edge_2id(pre_node, tmp_node, edge_iter);
        assert(find_edge);

        if (src == 15 and dest == 8)
        {
            //cout<<"find_edge: "<<find_edge<<endl;
            //cout<<"pre_node: "<<pre_node<<endl;
            //cout<<"target_node: "<<tmp_node<<endl;
            ////cout<<"Edge down_bw: "<< edge_iter->property().down_bw_residue <<endl;
            //cout<<"min_bw: "<<min_bw<<endl;
        }

        if (edge_iter->property().down_bw_residue < min_bw)
        {
            has_bw = 0;
            break;
        }

    }

    if (no_share_nodes && has_bw)
    {
        stack_nodes.insert(stack_nodes.end(),nodes_to_add.begin(),nodes_to_add.end());
        findPath = 1;
        return;
    }

    size_t cur_node = stack_nodes.back();
    size_t next_node;

    if (debug)
        cout<<"cur_code: "<<cur_node<<endl;

    vertex_iterator cur_vit = g.find_vertex(cur_node);
    cur_vit->property().visited = 1;

    vector<pair_IntFlt> & sorted_edges = cur_vit->property().sorted_edges_of_vertex;
    for(auto v_cost_pair = sorted_edges.begin(); v_cost_pair != sorted_edges.end(); ++v_cost_pair)
    {
        next_node = v_cost_pair->first;
        edge_iterator cur_edge_iter;

        //!!!!Can be optimized by modifying v_cost_pair to v_cost_bw_pair
        bool find_edge = g.find_out_edge_2id(cur_node, next_node, cur_edge_iter);
        assert(find_edge);

        has_bw = (cur_edge_iter->property().down_bw_residue >= min_bw);

        vertex_iterator next_vit = g.find_vertex(next_node);

        //push the next node only if the edge has enough bandwidth
        if ( (invalid_nodes.find(int(v_cost_pair->first)) == invalid_nodes.end()) &&  
             has_bw && !(next_vit->property().visited) )
        {
            stack_nodes.push_back(next_node);

            backTrackingSearch(g, invalid_nodes, src, dest, min_bw, stack_nodes, findPath);
            if (findPath)
                return;
            stack_nodes.pop_back();
        }
    }

    if (debug)
        cout<<"back track end."<<endl;
}
 


//find backup path
bool findOnePath(graph_t& g, size_t src, size_t dest, size_t min_bw, unordered_set<int>& invalid_nodes, 
    vector<size_t> &stack_nodes)
{
    size_t cur_node = src;
    bool find_path = 0;

    stack_nodes.push_back(cur_node);

    backTrackingSearch(g, invalid_nodes, src, dest, min_bw, stack_nodes, find_path);

    //cout<<"size of stack_nodes: "<<stack_nodes.size()<<endl;

    return find_path;
}


//shortest path with multi-constraints
void spmcSub(ofstream& myfile, graph_t& g, size_t src, size_t dest, size_t min_bw, size_t alpha \
            ,bool &find_path, vector<size_t> &shortest_path, vector<size_t> &backup_path, bool &recompute_g) // src and dest are exID
{
    //Run single src and dest
    if (recompute_g)
        sssp(g, src, dest, alpha);

    //cout << "!!!1"<<endl;

    //Compute reduced cost
    if (recompute_g)
        updateReducedCost(g, alpha);

    //cout << "!!!2"<<", src: "<<src<<", dest: "<<dest<<", min_bw :"<<min_bw<<endl;
    //Find shortest path
    //getShortestPath(g, src, dest, min_bw, shortest_path);
    unordered_set<int> empty_set = {-1};
    find_path |=  findOnePath(g, src, dest, min_bw, empty_set, shortest_path);

    //cout << "!!!3"<<", src: "<<src<<", dest: "<<dest<<", min_bw :"<<min_bw<<endl;

    //Find backup path
    unordered_set<int> invalid_nodes(shortest_path.begin(), shortest_path.end());
    //
    //cout<< "size of shortest path"<<shortest_path.size()<<", size of invalid_nodes "<<invalid_nodes.size()<<endl;
    find_path &= findOnePath(g, src, dest, min_bw, invalid_nodes, backup_path);

    //cout << "!!!4"<<endl;
}

//spmc: shortest path with multi-constraints
void spmc(ofstream& myfile, graph_t& g, vector<size_t> &src, vector<size_t> &min_bw, size_t dest,  \
          gBenchPerf_event & perf, int perf_group)
{

    perf.open(perf_group);
    perf.start(perf_group);
    myfile << "Src, Dst, Shortest Path, Back Up Path \n";

    #ifdef SIM
        SIM_BEGIN(true);  
    #endif

    //
    size_t max_alpha = pow(8,10);
    auto bw_iter = min_bw.begin();
    size_t alpha = 1;

    for (auto src_iter = src.begin(); src_iter != src.end(); ++ src_iter)
    {
        bool find_path = 0;
        bool reset_occurrence = 1, recompute_g;
        vector<size_t> shortest_path, backup_path;

        recompute_g = (src_iter == src.begin()) || (alpha != 1);
        alpha = 1;
    
        while (!find_path && alpha < max_alpha)
        {

            spmcSub(myfile, g, *src_iter, dest, *bw_iter, alpha, find_path, shortest_path, backup_path, recompute_g);

            //cout<<"Size of spa: "<<shortest_path.size()<<endl;
            //cout<<"Size of bpa: "<<backup_path.size()<<endl<<endl;

            alpha *= 40; //int(math.sqrt(rate));
            if (!find_path && alpha < max_alpha)
            {
                //cout<<"Value of find_path"<<endl;
                reset_occurrence = 0;
                recompute_g = 1;
                reset_graph(g, reset_occurrence);
                shortest_path.clear();
                backup_path.clear();
            }
        }

        assert(find_path);
        //if (find_path){
        //    //cout<<"Find the two shortest paths."<<endl;
        //    //cout<<"Size of shortest path: "<<shortest_path.size()<<endl;
        //    //cout<<"Size of backup path: "<<backup_path.size()<<endl<<endl;
        //}
        //else
        //    cout<<"Fail to find the two shortest paths.\n"<<endl;

        reset_occurrence = 1;
        reset_graph(g, reset_occurrence);
        advance(bw_iter, 1);

        //Print results
        myfile<<*src_iter<<","<<dest<<", ";
        for (auto iter = shortest_path.begin(); iter != shortest_path.end(); ++iter)
        {
            myfile<<*iter<<"|";
            //cout<<*iter<<"|";
        }

        myfile<<",";
        for (auto iter = backup_path.begin(); iter != backup_path.end(); ++iter)
        {
            myfile<<*iter<<"|";
            //cout<<*iter<<"|";
        }
    
        myfile<<"\n";
    }
    
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
    arg.add_arg("src","0","root/src vertex");
    arg.add_arg("dest","10","root/dest vertex");
}
//==============================================================//

class TEST_IN
{
    public:
        vector<vector<size_t> > src;
        vector<vector<size_t> > min_bw;
        vector<size_t> dest;

        void print()
        {
            auto iter_des = dest.begin();
            for (auto iter_i = src.begin(); iter_i != src.end(); ++iter_i)
            {
                cout<<"src nodes: ";
                for (auto iter_j = iter_i->begin(); iter_j != iter_i->end(); ++iter_j)
                {
                    cout<<*iter_j<<",";
                }
                cout<<" dest node: "<<*iter_des<<endl;
                advance(iter_des,1);
            }
        }
};

//Load AGG_pairs
TEST_IN loadTestIn(void)
{
    TEST_IN res;

   // new data
    size_t raw[21][3] = {{0,    1  , 10},
                      {0,    2  , 10},
                      {1  ,  2  , 10},
                      {1  ,  3  , 10},
                      {2  ,  3  , 10},
                      {0  ,  3  , 10},
                      {0  , 334 , 10},
                      {1  , 334 , 10},
                      {2  , 334 , 10},
                      {3  , 334 , 10},
                      {0  , 394 , 10},
                      {1  , 394 , 10},
                      {2  , 394 , 10},
                      {3  , 394 , 10},
                      {334, 394 , 10},
                      {0  , 557 , 10},
                      {2  , 557 , 10},
                      {1  , 557 , 10},
                      {3  , 557 , 10},
                      {334, 557 , 10},
                      {394, 557 , 10}};
    
    
    for (int i=0; i<21; i++)
    {
        vector<size_t> p_src;
        vector<size_t> p_min_bw;
        size_t p_dest = raw[i][1];

        if(!res.dest.empty() && res.dest.back() == p_dest)
        {
            res.src.back().push_back(raw[i][0]);
            res.min_bw.back().push_back(raw[i][2]);
        }
        else
        {
            p_src.push_back(raw[i][0]);
            p_min_bw.push_back(raw[i][2]);
            res.src.push_back(p_src);
            res.dest.push_back(p_dest);
            res.min_bw.push_back(p_min_bw);
        }
    }

    return res;
}


//Initialize a graph
void graphInit(graph_t &graph, string vfile, string efile, string separator)
{
    cout<<"loading data... \n";
    double t1 = timer::get_usec();
    
    if (graph.load_csv_vertices(vfile, true, ",", 0) == -1)
        return;
    if (graph.load_csv_edges(efile, true, ",", 1, 2,false, NULL) == -1) 
        return;

    size_t vertex_num = graph.num_vertices();
    size_t edge_num = graph.num_edges();
    double t2 = timer::get_usec();
    cout<<"== "<<vertex_num<<" vertices  "<<edge_num<<" edges\n";

    #ifndef ENABLE_VERIFY
        cout<<"== time: "<<t2-t1<<" sec\n\n";
    #endif

    // sanity check
    //uint64_t internel_srcID = graph.external_to_internel_id(to_string(src));
    //if (graph.find_vertex(internel_srcID)==graph.vertices_end()) 
    //
    //{
    //    cerr<<"wrong source vertex: "<<src<<endl;
    //    assert(false);
    //}
    //uint64_t internel_destID = graph.external_to_internel_id(to_string(dest));
    //if (graph.find_vertex(internel_destID)==graph.vertices_end()) 
    //{
    //    cerr<<"wrong dest vertex: "<<dest<<endl;
    //    assert(false);
    //}

}

//Generate several test pairs
TEST_IN genPairsIn(size_t total_vertex, size_t max_bw, size_t max_cost, size_t max_phyd)
{
    TEST_IN res;
    size_t num_pairs = 10000;

    vector<vector<size_t> > pairs;
    unordered_set<string> visited_pairs;

    for (size_t i=0; i<num_pairs; i++)
    {
        size_t p_src, p_dest, min_bw;
        vector<size_t> pair;
        string p_str;
        srand (1223);

        while (1)
        {
            p_src = rand() % (total_vertex);
            p_dest = rand() % (total_vertex);
            p_str = to_string(p_src) + "_" + to_string(p_dest);
            if (p_src != p_dest && (visited_pairs.find(p_str) == visited_pairs.end()) ) 
            {
                visited_pairs.insert(p_str);
                break;
            }
        }
        min_bw = 0;//rand() % (max_bw) + 1;
        pair.push_back(p_src); 
        pair.push_back(p_dest); 
        pair.push_back(min_bw); 
        pairs.push_back(pair);

    }

    //Get TEST_IN
    for (size_t i=0; i<pairs.size(); i++)
    {
        vector<size_t> p_src;
        vector<size_t> p_min_bw;
        size_t p_dest = pairs[i][1];

        if(!res.dest.empty() && res.dest.back() == p_dest)
        {
            res.src.back().push_back(pairs[i][0]);
            res.min_bw.back().push_back(pairs[i][2]);
        }
        else
        {
            p_src.push_back(pairs[i][0]);
            p_min_bw.push_back(pairs[i][2]);
            res.src.push_back(p_src);
            res.dest.push_back(p_dest);
            res.min_bw.push_back(p_min_bw);
        }
    }
    //Print res
    //res.print();

    return res;
}


double serialTest(string vfile, string efile, vector<vector<size_t> > &src, vector<vector<size_t> > &min_bw,  \
    vector<size_t> dest, string separator, gBenchPerf_event perf )
{
    graph_t graph;  
    graphInit(graph, vfile, efile, separator);

    double t1 = timer::get_usec();
    vector<size_t> p_src, p_min_bw;
    size_t p_dest;

    ofstream res_fstream;
    res_fstream.open("./test_res/serial.csv", ofstream::trunc);
     
    for (unsigned int i=0; i<src.size(); i++)
    {        
        p_src = src[i]; 
        p_dest = dest[i];
        p_min_bw = min_bw[i];

        //res_fstream.open("./test_res/parallel_res_"+to_string(idx), ofstream::trunc);
        spmc(res_fstream, graph, p_src , p_min_bw, p_dest, perf, i);  
    }

    res_fstream.close();
    double t2 = timer::get_usec();
    return t2-t1;
}



double parallelTest(string vfile, string efile, vector<vector<size_t> > &src, vector<vector<size_t> > &min_bw,  \
    vector<size_t> dest, string separator, gBenchPerf_event perf )
{
    //graph_t *graph;     
    double t1 = timer::get_usec();
    size_t thread_num = 64;
    assert(src.size() >= thread_num);

    //omp_set_num_threads(src.size());
    omp_set_num_threads(thread_num);
    #pragma omp parallel 
    {
        size_t i, idx, num_threads, num_tasks_per_thread;

        idx = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        //cout<<"number of threads: "<<num_threads<<endl;

        num_tasks_per_thread = src.size()/num_threads + (src.size()%num_threads > 0);

        ofstream res_fstream;
        res_fstream.open("./test_res/parallel_"+to_string(idx)+".csv");

        //graph = new graph_t();
        graph_t graph;
        graphInit(graph, vfile, efile, separator);

        for (i=idx*num_tasks_per_thread; i<(idx+1)*num_tasks_per_thread; i+=1)
        {
            if(i < src.size())
            {
                vector<size_t> p_src = src[i];    
                vector<size_t> p_min_bw = min_bw[i];
                size_t p_dest = dest[i];        
                spmc(res_fstream, graph, p_src , p_min_bw, p_dest, perf, 0); 
                //delete graph;
            }
        }
        double t_idx = timer::get_usec();
        res_fstream << "Run time: "<< t_idx-t1 <<endl;
        res_fstream.close();       
    }

    double t2 = timer::get_usec(); 
    return t2-t1;
}


int main(int argc, char * argv[])
{
    graphBIG::print();
    cout<<"Benchmark: shortest path with multiple constraints\n";

    argument_parser arg;
    gBenchPerf_event perf;
    arg_init(arg);
    if (arg.parse(argc,argv,perf,false)==false)
    {
        arg.help();
        return -1;
    }
    string path, separator;
    string vfile, efile;

    arg.get_value("dataset",path);
    arg.get_value("separator",separator);
     
    int total_vertex = 400, max_bw = 10, max_cost = 100, max_phyd = 100;
    vfile = path + "/vertex_"+to_string(total_vertex)+".csv";  
    efile = path + "/edge_"+to_string(total_vertex)+".csv";   
    
    cout<<"...\n";

    //Generate tests
    TEST_IN test = genPairsIn(total_vertex, max_bw, max_cost, max_phyd);
    int test_num = test.src.size();
    
    double s_time, p_time;

    cout<<"Total number of tests is "<< test_num <<endl;    

    //Serial run
    cout<<"Start running serial test"<<endl;
    s_time = serialTest(vfile, efile, test.src, test.min_bw, test.dest, separator, perf);

    #ifndef ENABLE_VERIFY
        cout<<"== Total ruuning time: "<<s_time<<" sec\n";
        cout<<"==================================================================\n"<<endl;
    #endif    

    //Parallel run
    cout<<"Start running openMP test"<<endl;
    p_time = parallelTest(vfile, efile, test.src, test.min_bw, test.dest, separator, perf);
    
    #ifndef ENABLE_VERIFY
        cout<<"== Total ruuning time: "<<p_time<<" sec\n";
        cout<<"==================================================================\n"<<endl;
    #endif    

    //Reseult verification
    #ifdef VERIFY_RESULTS
        bool all_identical = 1;
        for (unsigned int i=0; i<test_num; i++)
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

