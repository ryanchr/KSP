uint64_t  size_t  
graph_t   graph_Tau

// graph_t& g,  external id is thought as string when load csv.
uint64_t internel_src  = g.external_to_internel_id(to_string(src)); // external id --> internel id


//ori graph:
//vertex_property: min_cost,predecessor,sorted_edges_of_vertex(vector<pair_IntFlt>)
//edge_property:   cost, phy_dist, reduced_cost  

// from id to vertex_iterator, and versa vice.
vertex_iterator vit = g.find_vertex(internel_dest);   
size_t tmpId = vertex_iterator->id();
// from vertex_iterator to each vertex_property: min_cost,predecessor,sorted_edges_of_vertex(vector<pair_IntFlt>)
vit->property().min_cost = 0;  
// from vertex_iterator to edge_iterator (i.e., from vetex to each out edge)
for (edge_iterator eit = vit->edges_begin(); eit != vit->edges_end(); eit++)
    eit->property().cost = xx;  // from edge_iterator to each edge_property: cost, phy_dist, reduced_cost  
    size_t v = eit->target();                 // the vertex id of the out edge pointed to
    vertex_iterator v_vit = g.find_vertex(v); // the corresponding vertex_iterator, next can access the vertex property

//do not know vertex id or need to deal with each vertex, start by search each vertex
for (vertex_iterator vit=g.vertices_begin(); vit!=g.vertices_end(); vit++) 



//tau graph:
//vertex_property_Tau: at_KSPaths, predecessor, internel_id_of_g,    min_cost, sum_distance, sum_hops 
//edge_property_Tau: None

PQ_KSP_candidates_tau; // X, priority_queue with pair_FltInt, store the minCost and interenl ID of tau (the last id). 
KSPaths_lastID_tau // vector for output: store the top k shortest path id of tau.



class vertex_property // new
{
public:
    vertex_property():min_cost(FLT_MAX),predecessor(MY_INFINITY){}

    float min_cost;    // for shortest path
    uint64_t predecessor; // for shortest path
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


// tmp
// the first shortest path constructed in tau
    /*// construct the 1st shortest path in the pseudo-tree
    vertex_iterator src_vit_tau;
    vertex_iterator dest_vit_tau;
    edge_iterator eit_tau;
    vertex_iterator src_vit_g;
    vertex_iterator dest_vit_g;
    edge_iterator eit_g;

    src_vit_g = g.find_vertex(internel_src);  // internel_dest
    src_vit_tau = tau.add_vertex();
    src_vit_tau->property().at_KSPaths = true;//src_vit_tau->property().KSPaths_record.push_back(1); // this node is at the 1st shortest path.
    src_vit_tau->property().internel_id_of_g = src_vit_g->id(); // this is internel_src
    src_vit_tau->property().predecessor = MY_INFINITY;
    src_vit_tau->property().min_cost = 0; 
    src_vit_tau->property().sum_distance = 0; 
    src_vit_tau->property().sum_hops = 0;

    assert(internel_src == src_vit_g->id());
    size_t internel_src_tau = src_vit_tau->id();
    uint64_t tau_tmp_id = internel_src_tau;
    while (src_vit_g->id() != internel_dest)
    {
        dest_vit_g = g.find_vertex( src_vit_g->property().predecessor );
        dest_vit_tau = tau.add_vertex();
        find_result = g.find_out_edge_2id(src_vit_g->id(), dest_vit_g->id(), eit_g); // for eit_g->property().cost and eit_g->property().phy_dist
        assert(find_result);

        dest_vit_tau->property().at_KSPaths = true;//dest_vit_tau->property().KSPaths_record.push_back(1); // this node is at the 1st shortest path.
        dest_vit_tau->property().predecessor = tau_tmp_id; // 
        dest_vit_tau->property().internel_id_of_g = dest_vit_g->id();   
        dest_vit_tau->property().min_cost += eit_g->property().cost; //   
        dest_vit_tau->property().sum_distance += eit_g->property().phy_dist; //   
        dest_vit_tau->property().sum_hops += 1;

        tau_tmp_id = dest_vit_tau->id();  // for next iteration use
        tau.add_edge(src_vit_tau->id(),dest_vit_tau->id(),eit_tau); // note: put all info at vertex, see dest_vit_tau
        src_vit_g  = dest_vit_g;  // later check if right?
    } 
    PQ_KSP_candidates_tau.push( pair_FltInt(dest_vit_tau->property().min_cost, dest_vit_tau->id()) );
    KSPaths_lastID_tau.push_back(dest_vit_tau->id());*/



// the last part of MPS
        /*bool found_deviation_id = false; // the last node cannot be the deviation node
        while (tmp_id != internel_src_tau);
        {
            tmp_id = dest_vit_tau->property().predecessor;
            dest_vit_tau = tau.find_vertex(tmp_id);
            if (!found_deviation_id) 
            {
                pk_vkt_nodes.push(tmp_id);
            }
            //if (!found_deviation_id && !dest_vit_tau->property().KSPaths_record.empty())
            if (!found_deviation_id && dest_vit_tau->property().at_KSPaths==true)
            {
                deviation_id = tmp_id; // not used for this deviation_id. It is the top() of pk_vkt_nodes
                found_deviation_id = true;
                break;
            }
            dest_vit_tau->property().at_KSPaths = true; //dest_vit_tau->property().KSPaths_record.push_back(k); 
        } */
        //assert(tmp_id == internel_src_tau);
