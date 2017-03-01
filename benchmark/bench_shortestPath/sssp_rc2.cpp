
void top_ksp(graph_t& g, size_t src, size_t dest, size_t Kvalue, double max_phy_dist, size_t max_phy_hops, size_t min_phy_hops, gBenchPerf_event & perf, int perf_group, ofstream myfile)
{

    perf.open(perf_group);
    perf.start(perf_group);
#ifdef SIM
    SIM_BEGIN(true);  
#endif


    //ofstream myfile ("outputSingleLayer/combineMethods_single_0_2.csv"); // update the above number
    size_t trueMinCost_Iter = 5*1e6; // 1e6   (0-334:  1e7->152s)

    bool trueMinCost        = true;
    size_t curr_kValue      = 0;


    // first operate trueMinCost with trueMinCost_Iter iterations
    top_ksp_subFun(trueMinCost, trueMinCost_Iter, curr_kValue, myfile, g, src, dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops);

    // second, find the remaining (Kvalue - curr_kValue) paths by using reduced hops
    reset_graph(g);
    trueMinCost = false;
    top_ksp_subFun(trueMinCost, trueMinCost_Iter, curr_kValue, myfile, g, src, dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops);

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


double serialTest(string vfile, string efile, size_t src, size_t dest, string separator, int AGG_pair_num, int Kvalue, int max_phy_dist, int max_phy_hops, size_t min_phy_hops,  vector<vector<int> >tests, gBenchPerf_event perf)
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
		//res_fstream.open("./test_res/parallel_res_"+to_string(idx), ofstream::trunc);
        res_fstream.open("outputSingleLayer/combineMethods_single_0_2.csv");
		
        top_ksp(*graph, p_src , p_dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops, perf, 0, res_fstream);  
        delete graph;

        res_fstream.close();
    }

    double t2 = timer::get_usec();
    return t2-t1;
}



double parallelTest(string vfile, string efile, size_t src, size_t dest, string separator, int AGG_pair_num, int Kvalue, int max_phy_dist, int max_phy_hops, size_t min_phy_hops, vector<vector<int> >tests, gBenchPerf_event perf)
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
		res_fstream.open("outputSingleLayer/combineMethods_single_0_2.csv");

        graph = new graph_t();
        graphInit(*graph, vfile, efile, src, dest, separator);
		res_fstream.open("./test_res/parallel_res_"+to_string(idx), ofstream::trunc);

        top_ksp(*graph, p_src , p_dest, Kvalue, max_phy_dist, max_phy_hops, min_phy_hops, perf, 0, res_fstream);  
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
    arg.get_value("min_phy_hops",min_phy_hops);
    arg.get_value("threadnum",threadnum);

    
    AGG_pair_num = 24;   // 5    50    100 
    Kvalue = 30;         // 5   50    500  
    max_phy_dist = FLT_MAX;
    max_phy_hops = 9;
	min_phy_hops = 4;
	
    total_vertex = 400;  // 400   4000   10000 
    vfile = path + "/vertex400.csv";  // vertex400   vertex4000    vertex10000
    efile = path + "/vertex400dim8_edge.csv";       
	
	
    cout<<"...\n";
    gBenchPerf_multi perf_multi(threadnum, perf);

    //Generate tests
    vector<vector<int> > tests = genTests(total_vertex, AGG_pair_num);
    double s_time, p_time;

    cout<<"AGG_pair_num is "<< AGG_pair_num <<"; Kvalue is "<<Kvalue<<endl;
    cout<<"max_phy_dist is " << max_phy_dist <<endl;
    cout<<"max_phy_hops is " << max_phy_hops <<endl<<endl;
	cout<<"min_phy_hops is " << min_phy_hops <<endl;
    
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

