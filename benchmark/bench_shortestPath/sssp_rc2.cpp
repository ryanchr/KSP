
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
