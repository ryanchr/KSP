
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

