#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include <vector> 

static void parse(float &x, const char * s) {
    x = (float) atof(s);
}

// Removes \n from the end of line
inline void FIXLINE(char * s) {
    int len = (int) strlen(s)-1;
    if(s[len] == '\n') s[len] = 0;
}

struct edge_t {
    int from; 
    int to;
    float cost;
    float dist;

    int id;
    
    void add_edge (int _from, int _to, int _n, float _cost=0, float _dist=0) {
        from = _from;
        to = _to;
        cost = _cost;
        dist = _dist;
        id = _n;
    }
};


void convert_edgelist(std::string inputfile, bool has_weight = false) {
    
    if (has_weight) std::cout << "[LOG] Input file has weight!" << std::endl;

    FILE * inf = fopen(inputfile.c_str(), "r");

    size_t bytesread = 0;
    size_t linenum = 0;
    int nnz = 0;

    std::vector<edge_t> edges;

    if (inf == NULL) {
        std::cout << "[Error] Could not load :" << inputfile << " error: " << std::endl; 
    }
    assert(inf != NULL);
    
    std::cout << "[LOG] Reading in edge list format!" << std::endl;
   
    char s[1024];
    while(fgets(s, 1024, inf) != NULL) {
        linenum++;
        if (linenum % 10000000 == 0) {
            std::cout << "[LOG] Read " << linenum << " lines, " << bytesread / 1024 / 1024.  << " MB" << std::endl;
        }
        FIXLINE(s);
        bytesread += strlen(s);
        //if (s[0] == '#') continue; // Comment
        if (s[0] == ' ') continue; // Comment
        
        char delims[] = ",";
        char * t;
        float val;
        float val2;


        t = strtok(s, delims);
        if (t == NULL) {
            assert(false);
        }
        int id = atoi(t);

        t = strtok(NULL, delims);
        if (t == NULL) {
            assert(false);
        }
        int from = atoi(t);

        t = strtok(NULL, delims);
        if (t == NULL) {
            assert(false);
        }
        int to = atoi(t);
        nnz++;
        edge_t e; 
        
        if (has_weight) {
            t = strtok(NULL, delims);
            assert(t != NULL);
            parse(val, (const char*) t);
            
            t = strtok(NULL, delims);
            assert(t != NULL);
            parse(val2, (const char*) t);
            
            e.add_edge(from, to, id, val, val2);
        } else
            e.add_edge(from, to, id);
        
        edges.push_back(e);
 
    }

    fclose(inf);

    //////////std::string ss = inputfile + ".undirected";

    //std::string ss = "LogicalLink.csv";
    std::string ss = "PhysicalLink.csv";

    FILE * of = fopen(ss.c_str(), "w");
    if (has_weight) {
    	fprintf(of, " LinkID,SrcID, DstID,Cost,Distance\n");
        for(int i=0; i< edges.size(); i++) {
            edge_t e = edges[i];
            fprintf(of, "%d,%d,%d,%f,%f,\n", e.id, e.from, e.to, e.cost, e.dist);
            fprintf(of, "%d,%d,%d,%f,%f,\n", e.id, e.to, e.from, e.cost, e.dist);
        }
    
    } else { 
        fprintf(of, " LinkID,SrcID, DstID\n");
        for(int i=0; i< edges.size(); i++) {
            edge_t e = edges[i];
            fprintf(of, "%d,%d,%d,\n", e.id, e.from, e.to);
            fprintf(of, "%d,%d,%d,\n", e.id, e.to, e.from);
        }
    }

    fclose(of);
    std::cout << "[LOG] DONE!" << std::endl;

}

int main(int argc, char ** argv) {

    //std::string in = "LogicalLink_ori.csv"; // ./conversions.out
    std::string in = "PhysicalLink_ori.csv"; // ./conversions.out 1

    bool has_weight= ( argc > 1) ? true : false;
    convert_edgelist(in, has_weight);

    return 0;
}
