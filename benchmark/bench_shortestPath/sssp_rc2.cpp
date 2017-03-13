
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
