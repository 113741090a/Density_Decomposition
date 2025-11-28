#include <fstream>
#include <iostream>
#include <queue>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdio.h> 
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <unordered_map>
#include <unordered_set>
// #include <stdio.h>
// #include <stdlib.h>
#include <omp.h>

//to check memory usage
#include <sys/resource.h>
#include <malloc.h>

// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/boykov_kolmogorov_max_flow.hpp>
// using namespace boost;
//#include <lp_lib.h>
//#include </userhome/35/qxue/miniconda3/envs/mark/include/glpk.h>
//#include <glpk.h>

//#include <libglpk.h>
//using namespace std;
// #include <vector>
// #include <utility>

//#define NLINKS 500000000//25690705119//maximum number of edges for memory allocation, will increase if needed

#define NLINKS 5000000

typedef struct node node; // 提前声明结构体类型 node


struct node{
	unsigned id; //id of the node. 
  int old_id;
    double weight;
	  double r;//rate value of node id
    double r_correct; //correct value of the rate value;
    //double rho; //density of node id. r/weight.
    int *neighboredges;
    int degree;
    int part;
    //node *neighbors;
}; //One side of "node" represent hyperedge. The other side represent hyper-vertex. 

typedef struct {
	int s; //hyper edge.
    //node hyperedge;
	int t;
    //node hypernode;
  double a;//alpha value of edge (s,t), from s to t.
  //double b;
} link0;





//data structure we need to do the iterations. 
typedef struct {
	unsigned n;//number of hyper nodes
	unsigned long long m;//number of hyper edges
  unsigned long long e;//number of edges
	int *map;//correspondance between old and new nodeID
  int *map2;//correspondance between old and new nodeID, hyperedge
	int *newlabel;//correspondance between old and new nodeID
  int *newlabel2;//correspondance between old and new nodeID, hyperedge
  node *hyperedges;
  node *hypernodes;
	link0 *edges;//list of all edges, link.
	//node *nodes;//value associated to each node
	//node *nodes2;//value associated to each node
  //double *ne;//ne[i]=number of edges from i to nodes before (used for pava)
	//unsigned *cd;//cumulative degree
	//unsigned *cuts;
  unsigned iter;//number of iterations
  double tk;
} optim;

void Out_put_density(optim *opt, char *file_name);
void Out_put_weights(optim *opt, char *file_name);


//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

unsigned max2(unsigned a,unsigned b){
  a=(a>b) ? a : b;
	return a;
}


optim* readedgelist(char* edgelist){
	unsigned long long e1=NLINKS;
	optim *opt=(optim*)malloc(sizeof(optim));
	FILE *file;
	opt->n=0; //number of hypernodes
  opt->m=0; //number of hyperedges
	opt->e=0;
	file=fopen(edgelist,"r");
	//check if the file opened successfully
    if (file == NULL) {
        printf("Cannot open file: %s\n", edgelist);
        free(opt);
        exit(-1);
        //return NULL;
    }
	opt->edges=(link0*)malloc(e1*sizeof(link0) );
  unsigned s_tmp = 0;
  unsigned t_tmp = 0;
  //opt->edges[opt->e].s = s_tmp;
  //printf("No problem0 %llu \n", e1);
	while (fscanf(file,"%u %u", &(opt->edges[opt->e].s), &(opt->edges[opt->e].t))==2) {
    //opt->edges[opt->e].s = 
		opt->m=max2(opt->m,opt->edges[opt->e].s);
    //max3(opt->m,opt->edges[opt->e].s,1);
    opt->n=max2(opt->n,opt->edges[opt->e].t);
		if (opt->e++==e1) {
			e1+=NLINKS;
			opt->edges=(link0*)realloc(opt->edges,e1*sizeof(link0));
		}
    
	}
	fclose(file);
	opt->n++;
  opt->m++;
	opt->edges=(link0*)realloc(opt->edges,opt->e*sizeof(link0) );
	return opt;
}


//just relabel the vertices. Maybe not useful. 
void relabel(optim *opt) {
	int long long i,j, j2;
	// int *newlabel; //for hyper edge 

  // int *newlabel2; // for hyper node

  //hyperedge relabel.
	opt->newlabel=(int*)malloc(opt->m*sizeof(int));
	for (i=0;i<opt->m;i++) {
		opt->newlabel[i]=-1;
	}
	opt->map=(int*)malloc(opt->m*sizeof(int));
  //hypernode relabel
	opt->newlabel2=(int*)malloc(opt->n*sizeof(int));
	for (i=0;i<opt->n;i++) {
		opt->newlabel2[i]=-1;
	}
	opt->map2=(int*)malloc(opt->n*sizeof(int));

	j=0;
  j2 = 0;
  //printf("okay\n");
  //printf("%d, %d\n", opt->m, opt->n);
	for (i=0;i<opt->e;i++) {
    //printf("%d, %d, %d\n",i, j, j2);
		if (opt->newlabel[opt->edges[i].s]==-1){
			opt->newlabel[opt->edges[i].s]=j;
			opt->map[j]=opt->edges[i].s;
      opt->edges[i].s=opt->newlabel[opt->edges[i].s];
      j++;
		}
    else{
      opt->edges[i].s=opt->newlabel[opt->edges[i].s];
    }
		if (opt->newlabel2[opt->edges[i].t]==-1){
			opt->newlabel2[opt->edges[i].t]=j2;
			opt->map2[j2]=opt->edges[i].t;
      opt->edges[i].t=opt->newlabel2[opt->edges[i].t];
      j2++;
		}
    else{
      opt->edges[i].t=opt->newlabel2[opt->edges[i].t];
    }
	}
  //printf("okay3\n");
  opt->m = j;
	opt->n=j2;
	opt->map=(int*)realloc(opt->map,opt->m*sizeof(int));
  opt->map2=(int*)realloc(opt->map2,opt->n*sizeof(int));
}


//initialize the optim datastructure
void init(optim *opt){
	unsigned long long k;
  opt->iter=0;
	opt->hypernodes=(node*)malloc(opt->n*sizeof(node));
  opt->hyperedges=(node*)malloc(opt->m*sizeof(node));
  for (k=0;k<opt->n;k++){
    opt->hypernodes[k].id=k;
    opt->hypernodes[k].old_id=k;
    opt->hypernodes[k].r=0;
    opt->hypernodes[k].weight=1.0;
    opt->hypernodes[k].degree=0;
    //opt->hypernodes[k].neighboredges = (int *)malloc(opt->m*sizeof(int));
  }
  for (k=0;k<opt->m;k++){
    opt->hyperedges[k].id=k;
    opt->hyperedges[k].r=0;
    opt->hyperedges[k].weight=1.0;
    opt->hyperedges[k].degree=0; 
    //opt->hyperedges[k].neighboredges = (int *)malloc(opt->n*sizeof(int));
    //printf("idx = %d, original idx = %d  ", k, opt->map[k]);
  }
	for (k=0;k<opt->e;k++){
    opt->hyperedges[opt->edges[k].s].degree++;
    opt->hypernodes[opt->edges[k].t].degree++;
	}
  for (k=0;k<opt->n;k++){
    opt->hypernodes[k].neighboredges = (int *)malloc(opt->hypernodes[k].degree*sizeof(int));
    opt->hypernodes[k].degree=0;
  }
  for (k=0;k<opt->m;k++){
    opt->hyperedges[k].neighboredges = (int *)malloc(opt->hyperedges[k].degree*sizeof(int));
    opt->hyperedges[k].degree=0;
  }
  //------read edge part-----------------------
  // FILE *file;
  
	// file=fopen("bipartite_mark_E.txt","r");
  // int idx;
  // double weight;
	// while (fscanf(file,"%d %lf", &(idx), &(weight))==2) {
  //   //opt->edges[opt->e].s = 
	// 	int new_idx = opt->newlabel[idx];
  //   opt->hyperedges[new_idx].weight = weight;
  // }
	// fclose(file);

	// file=fopen("bipartite_mark_N.txt","r");
	// while (fscanf(file,"%d %lf", &(idx), &(weight))==2) {
  //   //opt->edges[opt->e].s = 
	// 	int new_idx = opt->newlabel2[idx];
  //   opt->hypernodes[new_idx].weight = weight;
  // }
	// fclose(file);



  //---------------------------------------
  //#pragma omp parallel for private(k)
	for (k=0;k<opt->e;k++){
    int hyperedge_idx = opt->edges[k].s;
    double weight = opt->hyperedges[hyperedge_idx].weight;
    opt->edges[k].a= weight / 2.0;
    //opt->edges[k].a= 1.0 / ((double) opt->n); //weighted, set to 1/2 for FISTA
		//#pragma omp atomic update
		opt->hyperedges[opt->edges[k].s].r= opt->hyperedges[opt->edges[k].s].r + opt->edges[k].a;
    
    opt->hyperedges[opt->edges[k].s].neighboredges[opt->hyperedges[opt->edges[k].s].degree] = k;
    opt->hyperedges[opt->edges[k].s].degree++;
    //#pragma omp atomic update
		opt->hypernodes[opt->edges[k].t].r= opt->hypernodes[opt->edges[k].t].r + opt->edges[k].a;
    //printf("access hypernode\n");
    opt->hypernodes[opt->edges[k].t].neighboredges[opt->hypernodes[opt->edges[k].t].degree] = k;
    opt->hypernodes[opt->edges[k].t].degree++;
	}
  // for (k=0;k<opt->n;k++){
  //   opt->hypernodes[k].neighboredges = (int *)realloc(opt->hypernodes[k].neighboredges, opt->hypernodes[k].degree*sizeof(int));
  // }
  // for (k=0;k<opt->m;k++){
  //   opt->hyperedges[k].neighboredges = (int *)realloc(opt->hyperedges[k].neighboredges, opt->hyperedges[k].degree*sizeof(int));
  // }
  //--------------------maybe delete later
  // int i = 0;
  // for (i=0; i < opt->n; i++){
  //   opt->hypernodes[i].r = 0;
  // }
  // for (i=0; i < opt->e; i++){
  //   int hyper_node_idx = opt->edges[i].t;
  //   opt->hypernodes[hyper_node_idx].r = opt->hypernodes[hyper_node_idx].r + opt->edges[i].a;
  // }
  // Out_put_density(opt, "output_init2.txt");
  // Out_put_weights(opt, "output_init2w.txt");
}
void update_density2(optim *opt, link0* weight);

//proportional projection
void onepass_PR_plus_product_normal_2(optim *opt, link0* last_x, link0* z){
	int i,j;
	unsigned long long k;

	opt->iter++;
  int t = opt->iter;
  //copy current t to last t
  for (k = 0; k < opt->e; k++){
    last_x[k].a = opt->edges[k].a;
  }

  double payload;
	//do PR on z. 
  update_density2(opt, z);

  //--do PR on z. 
  for (i=0;i<opt->m;i++){
    payload = 0.0;
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double r = opt->hypernodes[neighbor_node_idx].r;
      
      payload = payload + opt->hypernodes[neighbor_node_idx].weight * z[edge_idx].a / r;
    }
//be carefull here. 
    if (payload == 0.0){
      for (j=0; j < opt->hyperedges[i].degree; j++){
        int edge_idx = opt->hyperedges[i].neighboredges[j];
        int neighbor_node_idx = opt->edges[edge_idx].t;
        double r = opt->hypernodes[neighbor_node_idx].r;
        opt->edges[edge_idx].a = opt->hyperedges[i].weight / opt->hyperedges[i].degree;
      }
      continue;
    }
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double r = opt->hypernodes[neighbor_node_idx].r;
      opt->edges[edge_idx].a = opt->hyperedges[i].weight * (opt->hypernodes[neighbor_node_idx].weight * z[edge_idx].a / r)/payload;
    }
  }

  //linear combination. 
  //double lambda1 = (1 + sqrt(1 + 4 * pow(*lambda, 2))) /2.0;
  //double gamma = (double)(t-1)/(double)(t+2);
  double gamma = (double)(t-1)/(double)(t+2);
  for (k = 0; k < opt->e; k++){
    //printf("values %lf, %lf\n", opt->edges[k].a, last_x[k].a);
    if (opt->edges[k].a == 0.0 || last_x[k].a == 0.0) z[k].a = 0.0;
    else z[k].a = pow(opt->edges[k].a, 1+gamma) * pow(last_x[k].a, -gamma);
    //exp((1+gamma)*log(opt->edges[k].a) - gamma*log(last_x[k].a));
  }


  //do projection. proportionally
  //#pragma omp parallel for private(i,j,k)
  for (i=0; i < opt->m; i++){
    double hyperedge_weight = opt->hyperedges[i].weight;
    double total_weight = 0.0;
    for (int i1 = 0; i1 < opt->hyperedges[i].degree; i1++){
      int edge_idx = opt->hyperedges[i].neighboredges[i1];
      total_weight = total_weight + z[edge_idx].a;
    }
    double parameter = hyperedge_weight / (total_weight);
    for (int i1 = 0; i1 < opt->hyperedges[i].degree; i1++){
      int edge_idx = opt->hyperedges[i].neighboredges[i1];
      z[edge_idx].a = z[edge_idx].a * parameter;
    }
  }
}


//one pass over all edges
void onepass_PR(optim *opt){
	unsigned i,j;
	unsigned long long k;
	double gamma;
	opt->iter++;
  double payload;
	//#pragma omp parallel for private(i,j,k)
  for (i=0;i<opt->m;i++){
    payload = 0.0;
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double r = opt->hypernodes[neighbor_node_idx].r;
      
      payload = payload + opt->hypernodes[neighbor_node_idx].weight * opt->edges[edge_idx].a / r;
    }
//be carefull here. 
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double r = opt->hypernodes[neighbor_node_idx].r;
      opt->edges[edge_idx].a = opt->hyperedges[i].weight * (opt->hypernodes[neighbor_node_idx].weight * opt->edges[edge_idx].a / r)/payload;
    }
  }
  for (i=0; i < opt->n; i++){
    opt->hypernodes[i].r = 0.0;
  }
  for (i=0; i < opt->e; i++){
    int hyper_node_idx = opt->edges[i].t;
    opt->hypernodes[hyper_node_idx].r = opt->hypernodes[hyper_node_idx].r + opt->edges[i].a;
  }
  //printf("one iteration\n");
  // for (i=0;i<opt->m;i++){
  //   for (j=0; j < opt->hyperedges[i].degree; j++){
  //     int edge_idx = opt->hyperedges[i].neighboredges[j];
  //     int neighbor_node_idx = opt->edges[edge_idx].t;
  //     opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + opt->edges[edge_idx].a;
  //   }
  // }
}



void Out_put_weights(optim *opt, char *file_name){
  FILE *file;
  file = fopen(file_name, "w"); // 以写入模式打开文件
  printf("%s", file_name);
  if (file == NULL) {
      printf("Cannot Open\n");
      return;
  }
  int i = 0;
  for (i=0; i < opt->e; i++){
    fprintf(file, "%d, %f\n", i, opt->edges[i].a);
  }
  fclose(file); 
  
}

void Out_put_density(optim *opt, char *file_name){
  FILE *file;
  file = fopen(file_name, "w"); // 以写入模式打开文件
  printf("%s", file_name);
  if (file == NULL) {
      printf("Cannot Open\n");
      return;
  }
  int i = 0;
  for (i=0; i < opt->n; i++){
    fprintf(file, "%d, %f\n", i, opt->hypernodes[i].r);
  }
  fclose(file); 
  
}

//used for quicksort, smaller r has higher priority. Be careful here.  
static int compare_nodes(void const *a, void const *b){
	node const *pa = (const node*)a;
	node const *pb = (const node*)b;
	if ((*pa).r<(*pb).r)
		return 1;
  else if ((*pa).r>(*pb).r){
    return -1;
  }
  else if ((*pa).id<(*pb).id){
    return -1;
  }
  else return 1;
	return -1;
}


int hasIntersection(int arr1[], int size1, int arr2[], int size2) {
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            if (arr1[i] == arr2[j]) {
                return 1; // 有交集
            }
        }
    }
    return 0; // 没有交集
}
//may be changed to weighted version later. 
int get_max_degree(optim *opt){
  int max_degree = 0;
  for (int i=0; i < opt->n; i++){
    if(opt->hypernodes[i].degree > max_degree){
      max_degree = opt->hypernodes[i].degree;
    } 
   // opt->hypernodes[i].degree = 0;
  }
  return max_degree;
}

double get_max_weighted_degree(optim *opt){
  double max_degree = 0.0;
  for (int i=0; i < opt->n; i++){
    if(opt->hypernodes[i].degree / opt->hypernodes[i].weight > max_degree){
      max_degree = opt->hypernodes[i].degree / opt->hypernodes[i].weight;
    } 
   // opt->hypernodes[i].degree = 0;
  }
  return max_degree;
}

//update density vector
void update_density(optim *opt){
  for (int i=0; i < opt->n; i++){
    opt->hypernodes[i].r = 0.0;
  }
  for (int i=0; i < opt->e; i++){
    int hyper_node_idx = opt->edges[i].t;
    opt->hypernodes[hyper_node_idx].r = opt->hypernodes[hyper_node_idx].r + opt->edges[i].a;
  }
}

//update density vector 2
void update_density2(optim *opt, link0* weight){
  for (int i=0; i < opt->n; i++){
    opt->hypernodes[i].r = 0.0;
  }
  for (int i=0; i < opt->e; i++){
    int hyper_node_idx = weight[i].t;
    opt->hypernodes[hyper_node_idx].r = opt->hypernodes[hyper_node_idx].r + weight[i].a;
  }
}


using namespace std;

#include <boost/config.hpp>
 
 
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
 
 
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>
 
 
/*
,
boost::property <boost::vertex_color_t, boost::default_color_type,
boost::property <boost::vertex_distance_t, long,
boost::property <boost::vertex_predecessor_t, Traits::edge_descriptor > > >
*/
typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
	boost::property < boost::vertex_name_t, std::pair<int,int>,
	boost::property < boost::vertex_color_t, boost::default_color_type,
	boost::property < boost::vertex_distance_t, long,
	boost::property < boost::vertex_predecessor_t, Traits::edge_descriptor >  > > >,
	boost::property < boost::edge_capacity_t, double,
	boost::property < boost::edge_residual_capacity_t, double,
	boost::property < boost::edge_reverse_t, Traits::edge_descriptor> > >
> Graph;
 
 
void AddEdge(	Traits::vertex_descriptor &v1,
				Traits::vertex_descriptor &v2,
				boost::property_map < Graph, boost::edge_reverse_t >::type &rev,
				const double capacity, 
				Graph &g)
{
	Traits::edge_descriptor e1 = boost::add_edge(v1, v2, g).first;
	Traits::edge_descriptor e2 = boost::add_edge(v2, v1, g).first;
	boost::put(boost::edge_capacity, g, e1, capacity);
	boost::put(boost::edge_capacity, g, e2, 0);
 
	rev[e1] = e2;
	rev[e2] = e1;
	
}
 


//verify use maximum flow
typedef struct {
	int n;//number of hypernodes
	int m;//number of hyperedges
	unsigned e;//number of links(edges)
  node *hyperedges;
  node *hypernodes;
	link0 *edges;//list of all edges, link.
} subgraph;

subgraph *allocsubgraph(int n, int m, unsigned e){
	subgraph *sg=(subgraph*)malloc(sizeof(subgraph));
  //printf("okay\n");
	sg->n=n;
  sg->m = m;
	sg->e=e;
	//sg->hyperedges=(node*)calloc(m,sizeof(node));
  sg->hyperedges=(node*)malloc(m*sizeof(node));
  //printf("okay\n");
	//sg->hypernodes=(node*)calloc(n,sizeof(node));
  sg->hypernodes=(node*)malloc(n*sizeof(node));
  //printf("okay\n");
	sg->edges=(link0*)malloc(e*sizeof(link0));
  //printf("Subgraph initialized\n");
	return sg;
}

void freesubgraph(subgraph* sg){
  // for (int i = 0; i < sg->m; i++){
  //   if (sg->hyperedges[i].neighboredges != NULL) {
  //     free(sg->hyperedges[i].neighboredges);      // 释放内存
  //     sg->hyperedges[i].neighboredges = NULL;     // 将指针设为 NULL，防止悬空指针
  //   }
  // }
	free(sg->hyperedges);
  // for (int i = 0; i < sg->n; i++){
  //   if (sg->hypernodes[i].neighboredges != NULL) {
  //     free(sg->hypernodes[i].neighboredges);      // 释放内存
  //     sg->hypernodes[i].neighboredges = NULL;     // 将指针设为 NULL，防止悬空指针
  //   }
  // }
	free(sg->hypernodes);
	free(sg->edges);
	free(sg);
}

void freesubgraph2(subgraph* sg){
  for (int i = 0; i < sg->m; i++){
    if (sg->hyperedges[i].neighboredges != NULL) {
      free(sg->hyperedges[i].neighboredges);      // 释放内存
      sg->hyperedges[i].neighboredges = NULL;     // 将指针设为 NULL，防止悬空指针
    }
  }
	free(sg->hyperedges);
  for (int i = 0; i < sg->n; i++){
    if (sg->hypernodes[i].neighboredges != NULL) {
      free(sg->hypernodes[i].neighboredges);      // 释放内存
      sg->hypernodes[i].neighboredges = NULL;     // 将指针设为 NULL，防止悬空指针
    }
  }
	free(sg->hypernodes);
	//free(sg->edges);
	free(sg);
}


static double threshold = 0.000000001;

void maxflow(subgraph* sg,FILE *file, node* Mark_nodes) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  long memory_usage_mb = usage.ru_maxrss / 1024;
  printf("Memory usage: %ld MB\n", memory_usage_mb);
  double hyperedge_weights = 0.0;
  double hypernode_weights = 0.0;
  for (int i = 0; i < sg->m; i++){
    hyperedge_weights = hyperedge_weights + sg->hyperedges[i].weight;
  }
  int* idx_map = (int*)malloc(sg->n*sizeof(int));
  for (int i = 0; i < sg->n; i++){
    hypernode_weights = hypernode_weights + sg->hypernodes[i].weight;
    idx_map[i] = -1;
  }
  double parameter = hyperedge_weights / hypernode_weights;
  //long int V = num_hyperedges+num_hypernodes+2; // Number of vertices
  Graph g;
	boost::property_map < Graph, boost::edge_reverse_t >::type rev = get(boost::edge_reverse, g);
  
  //printf("initialize graph sucessfull\n");
  //int s = 0;
  Graph::vertex_descriptor s = boost::add_vertex(g);
  Graph::vertex_descriptor t = boost::add_vertex(g);
  int hyper_idx = 2;
  for (int k1 = 0; k1 < sg->m; k1++){
    Graph::vertex_descriptor u = boost::add_vertex(g);
    hyper_idx++;
    AddEdge(s, u, rev, sg->hyperedges[k1].weight, g);
    //g.addEdge(s, idx, opt->hyperedges[k1].weight);
    for (int k2 = 0; k2 < sg->hyperedges[k1].degree; k2++){
      int edge_idx = sg->hyperedges[k1].neighboredges[k2];
      int original_hyper_node_idx = sg->edges[edge_idx].t;
      Graph::vertex_descriptor v;
      if (idx_map[original_hyper_node_idx]<0){
        v= boost::add_vertex(g);
        idx_map[original_hyper_node_idx] = hyper_idx++;
        AddEdge(v, t,rev, parameter*sg->hypernodes[original_hyper_node_idx].weight, g);
      }
      else{
        v = boost::vertex(idx_map[original_hyper_node_idx],g);
      }
      //g.addEdge(idx, new_idx, opt->hyperedges[k1].weight);
      AddEdge(u, v,rev, sg->hyperedges[k1].weight, g);
    }
  }
  
  // int t = num_hyperedges + num_hypernodes + 1;
  // for (int k1 = 0; k1 + i1 <= i2; k1++){
  //   int node_idx = node_arr[k1+i1].id;
  //   //g.addEdge(k1 + num_hyperedges + 1, t, parameter*node_arr[k1+i1].weight);
  //   add_edge(k1 + num_hyperedges + 1, t, parameter*node_arr[k1+i1].weight, g);
  // }
//  printf("add all edges graph sucessfull\n");
  //double flow = push_relabel_max_flow(g, s, t);
	double flow = boost::boykov_kolmogorov_max_flow(g, s, t);
	//double flow = boost::push_relabel_max_flow(g, s, t);
	//double flow = boost::edmonds_karp_max_flow(g, s, t);
//  printf("hyeredgeweights = %lf, flow value = %lf \n", hyperedge_weights, flow);
  
  //double threshold = 0.000000001;
  double z = flow;
  boost::property_map < Graph, boost::edge_residual_capacity_t >::type mark = get(boost::edge_residual_capacity, g);
	boost::property_map < Graph, boost::vertex_color_t >::type col = get(boost::vertex_color, g);


  // 遍历所有边并访问它们 boost::graph_traits < Graph >::
  // boost::graph_traits < Graph >::edge_iterator ei, ei_end;
  // for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
  //     Traits::edge_descriptor current_edge = *ei;
  //     Graph::vertex_descriptor u = source(current_edge, g);
  //     Graph::vertex_descriptor v = target(current_edge, g);
  //     if (v==t && mark[current_edge] == 0) col[u] = boost::red_color;
  //     std::cout << "Edge (" << source(current_edge, g) << " -> " << target(current_edge, g) << ")" << std::endl;
  //     std::cout << "Edge (" << mark[current_edge] << std::endl;
  // }


  // if (z -  (hyperedge_weights) < threshold && z - (hyperedge_weights) > - threshold){
  //   fprintf(file,"%u %le\n",sg->n,((double)(parameter)));
  //   return;
  // }
  if (z -  (hyperedge_weights) < threshold && z - (hyperedge_weights) > - threshold){
  //  printf("separated\n");
    fprintf(file,"%u %le\n",sg->n,((double)(parameter)));
    for (int x = 0; x < sg->n; x++){
      int idx = sg->hypernodes[x].old_id;
      //printf("idx = %d is modified to r_correct = %f\n", idx, parameter);
      Mark_nodes[idx].r_correct = parameter;
      //printf("Mark_nodes[idx].r_correct = %f\n", Mark_nodes[idx].r_correct);
    }
    //----------------------------------
    // FILE *file_density;
    // file_density = fopen("file_density.txt", "a"); // 以写入模式打开文件
    // if (file_density == NULL) {
    //   printf("Cannot Open\n");
    // }
    // int i = 0;
    // for (int x = 0; x < sg->n; x++){
    //   int idx = sg->hypernodes[x].old_id;
    //   fprintf(file_density, "%d, %.16f\n", idx, parameter);
    //   //printf("Mark_nodes[idx].r_correct = %f\n", Mark_nodes[idx].r_correct);
    // }
    // fclose(file_density); 
    //----------------------
    g.clear();
    freesubgraph(sg);
    free(idx_map);
    return;
  }

	bool *isin=(bool *)calloc(sg->n,sizeof(bool));
	unsigned n_nodes=0;
  int graph_idx = 0;
	for(unsigned i=0;i<sg->n;i++){
    graph_idx = idx_map[i];
		if (col[vertex(graph_idx,g)]==4){ //full
			isin[i]=1;
			n_nodes++;
			
		}
    //std::cout << i << " " << col[vertex(graph_idx,g)] << std::endl;
	}

	g.clear();
  //----min-cut suggest that the part is full but not maxflow
  if (n_nodes == sg->n || n_nodes == 0){
  //  printf("flow not full but do not separate, use mincut seperate, diff = %lf\n", z - (hyperedge_weights));
    fprintf(file,"%u %le\n",sg->n,((double)(parameter)));
    for (int x = 0; x < sg->n; x++){
      int idx = sg->hypernodes[x].old_id;
      //printf("idx = %d is modified to r_correct = %f\n", idx, parameter);
      Mark_nodes[idx].r_correct = parameter;
      //printf("Mark_nodes[idx].r_correct = %f\n", Mark_nodes[idx].r_correct);
    }
    //----------------------------------
    // FILE *file_density;
    // file_density = fopen("file_density.txt", "a"); // 以写入模式打开文件
    // if (file_density == NULL) {
    //   printf("Cannot Open\n");
    // }
    // int i = 0;
    // for (int x = 0; x < sg->n; x++){
    //   int idx = sg->hypernodes[x].old_id;
    //   fprintf(file_density, "%d, %.16f\n", idx, parameter);
    //   //printf("Mark_nodes[idx].r_correct = %f\n", Mark_nodes[idx].r_correct);
    // }
    // fclose(file_density); 
    //----------------------
    freesubgraph(sg);
    free(idx_map);
    free(isin);
    return;
  }

  subgraph *sg1=allocsubgraph(n_nodes,sg->m, sg->e);
  subgraph *sg2=allocsubgraph(sg->n-n_nodes,sg->m, sg->e);
  unsigned *newlabel=(unsigned*)malloc(sg->n*sizeof(unsigned));
  unsigned j1=0,j2=0;
  for(unsigned i=0;i<sg->n;i++){
    if (isin[i]){
      sg1->hypernodes[j1] = sg->hypernodes[i];
      newlabel[i]=j1++;
    }
    else{
      sg2->hypernodes[j2] = sg->hypernodes[i];
      newlabel[i]=j2++;
    }
  }
  sg1->n = j1;
  sg2->n = j2;
  int j1_edge = 0;
  int j2_edge = 0;
  j1=0, j2=0;
  for (int i = 0; i < sg->m; i++){
    int flag = 0;
    for (int j = 0; j < sg->hyperedges[i].degree; j++){
      int edge_idx = sg->hyperedges[i].neighboredges[j];
      if (! isin[sg->edges[edge_idx].t]){
        flag = 1;
      }
    }
    if (flag == 0){
      sg1->hyperedges[j1] = sg->hyperedges[i];
      for (int j = 0; j < sg->hyperedges[i].degree; j++){
        int edge_idx = sg->hyperedges[i].neighboredges[j];
        sg1->edges[j1_edge].s = j1;
        sg1->edges[j1_edge].t = newlabel[sg->edges[edge_idx].t];
        sg1->hyperedges[j1].neighboredges[j] = j1_edge;
        j1_edge++;
      }
      j1++;
    }
    else{
      sg2->hyperedges[j2] = sg->hyperedges[i];
      int j_tmp = 0;
      for (int j = 0; j < sg->hyperedges[i].degree; j++){
        int edge_idx = sg->hyperedges[i].neighboredges[j];
        if (! isin[sg->edges[edge_idx].t]){
          sg2->edges[j2_edge].s = j2;
          sg2->edges[j2_edge].t = newlabel[sg->edges[edge_idx].t];
          sg2->hyperedges[j2].neighboredges[j_tmp] = j2_edge;
          j2_edge++;
          j_tmp++;
        }
      }
      sg2->hyperedges[j2].degree = j_tmp; 
      j2++;
    }
  }
  sg1->m = j1;
  sg2->m = j2;
  sg1->e=j1_edge;
  sg2->e=j2_edge;
  free(isin);
  freesubgraph(sg);
  free(newlabel);
  free(idx_map);
  maxflow(sg1,file, Mark_nodes);
  maxflow(sg2,file, Mark_nodes);

}



int verify_maxflow(subgraph* sg, node* Mark_nodes){
  FILE *file;
  file = fopen("maxflow_res_large.txt", "a"); // 以写入模式打开文件
  printf("Start to verify using maxflow\n");
  maxflow(sg, file, Mark_nodes);
  fclose(file);
  return 1;
}

void freeopt(optim *opt){
  free(opt->map);
  free(opt->map2);
  free(opt->newlabel);
  free(opt->newlabel2);
  for (int i = 0; i < opt->m; i++){
    if (opt->hyperedges[i].neighboredges != NULL) {
      free(opt->hyperedges[i].neighboredges);      // 释放内存
      opt->hyperedges[i].neighboredges = NULL;     // 将指针设为 NULL，防止悬空指针
    }
  }
  free(opt->hyperedges);
  for (int i = 0; i < opt->n; i++){
    if (opt->hypernodes[i].neighboredges != NULL) {
      free(opt->hypernodes[i].neighboredges);      // 释放内存
      opt->hypernodes[i].neighboredges = NULL;     // 将指针设为 NULL，防止悬空指针
    }
  }
  free(opt->hypernodes);
  free(opt->edges);
}

void read_edgeweights(optim *opt, char *name1){
  //------read edge part-----------------------
  FILE *file;
  
	file=fopen(name1,"r");
  unsigned idx;
  double weight;
  // 使用 fscanf 跳过一行
  fscanf(file, "%u", &idx);
  while (fscanf(file,"%u, %lf", &(idx), &(weight))==2) {
    //opt->edges[opt->e].s = 
    opt->edges[idx].a = weight;
  }
	fclose(file);

}

int verify(optim* opt, node* Mark_nodes){
  //used for sorting
  node* node_arr = (node*)malloc(opt->n * sizeof(node));
  int count = 0;

  for (int i = 0; i < opt->n; i++){
    node_arr[i].id = opt->hypernodes[i].id;
    node_arr[i].r = opt->hypernodes[i].r;
    node_arr[i].degree = opt->hypernodes[i].degree;
    node_arr[i].weight = opt->hypernodes[i].weight;
  }
  qsort(node_arr,opt->n,sizeof(node),compare_nodes);
  //int* map = (int*)malloc(opt->n * sizeof(int)); //remember release this.
  //printf("%f ", node_arr[i].r);
  for (int i = 0; i<opt->n; i++){
    //printf("node id is %d \n", node_arr[i].id);
    opt->map2[node_arr[i].id] = i;
    //map[node_arr[i].id] = i;
    opt->hypernodes[i].r = 0;
  }
  //now abuse field r. r is y. 
  //every hyperedge gives its weight to its neighbor with largest label
  for (int i = 0; i < opt->m; i++){
    int max_idx = 0;
    int j_idx = 0;
    for (int j = 0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int original_idx = opt->edges[edge_idx].t;
      int new_idx = opt->map2[original_idx];
      if (new_idx >= max_idx){
        max_idx = new_idx;
        j_idx = j;
      }
    }
    int edge_idx = opt->hyperedges[i].neighboredges[j_idx];
    int original_idx = opt->edges[edge_idx].t;
    opt->hypernodes[original_idx].r = opt->hypernodes[original_idx].r + opt->hyperedges[i].weight; //weight... change to weighted version later.
  }
  double* y_arr = (double*)malloc(opt->n*sizeof(double));
  for (int i = 0; i<opt->n; i++){
    int idx = node_arr[i].id;
    y_arr[i] = opt->hypernodes[idx].r;
  }
  //Now use pava algorithm
	int *nag=(int*)malloc(opt->n*sizeof(int)); //number of elements of a part
	double *val=(double*)malloc(opt->n*sizeof(double)); //average value of a part
	int i,j;

	nag[0]=1;
	val[0]=y_arr[0];
	j=0;
	for (i=1;i<opt->n;i++){
		j+=1;
		val[j]=y_arr[i];
		nag[j]=1;
		while ((j>0) && (val[j]>val[j-1]-1e-10)){//do val[j]>val[j-1] to have a non-increasing monotonic regression.
			val[j-1]=(nag[j]*val[j]+nag[j-1]*val[j-1])/(nag[j]+nag[j-1]);
			nag[j-1]+=nag[j];
			j--;
		}
	}
  int num_componenets = j+1;
  //tentative decomposition
  int *decom=(int*)malloc(num_componenets*sizeof(int));
  int num = 0;
  //Now try to decompose. 
  int sep_idx = 0;
  int acc = 0; //accumulated idx.
  printf("Print parts\n");
  for (int k = 0; k < num_componenets; k++){
    //sep_idx = sep_idx + nag[k] - 1;
    acc = acc + nag[k];
    sep_idx = acc-1;
    decom[num++] = acc;
    //fprintf(file_arr, "%d\n", acc);
  }
  //fclose(file_arr);
  free(y_arr);
  free(val);



  //check stability
  printf("Check stability, %d\n", num_componenets);
  update_density(opt);

  bool *mark_hyper = (bool*)malloc(opt->m*sizeof(bool));
  std::vector<node> nodes;
  double* temp_r = (double*)malloc(opt->n*sizeof(double));
  for (int k =0; k < opt->n; k++){
    temp_r[k] = opt->hypernodes[k].r;
  }
  for (int k =0; k < opt->m; k++){
    mark_hyper[k] = 0;
    //nodeMap[opt->hyperedges[k].id] = opt->hyperedges[k];
  }
  int mark;
  //scanf("slkdjflkdjf\n %d", &mark);
  acc = 0;
  sep_idx = 0;
  num = 0;
  //FILE *file_arr = fopen("file_arr4.txt", "w");
  std::unordered_set<int> mySet; //store hyperedges contains crossing links
  for (int k = 0; k < num_componenets; k++){
    printf("part number is %d\n", nag[k]);
    int i1 = acc;
    acc = acc + nag[k];
    sep_idx = acc-1;
    int i2 = sep_idx;
    //printf("%d, %d", acc, sep_idx);
    //printf("Still works here, %d\n", num_componenets);
    int flag_k = 0;
    for (int i = i1; i <= i2; i++){
      int node_idx = node_arr[i].id;
      //if (i == 416) printf("node_idx = %d\n", node_idx);
      for (int j=0; j < opt->hypernodes[node_idx].degree; j++){
        int edge_idx = opt->hypernodes[node_idx].neighboredges[j];
        int hyperedge_idx = opt->edges[edge_idx].s;
        
        int num_edges = 0;
        int degree = opt->hyperedges[hyperedge_idx].degree;
        //int *mark_edges = (int*)malloc(degree*sizeof(int));
        double ratio = 0.0;
        double total_payload = 0.0;
        double mark_payload = 0.0;
        sep_idx = i2;
        for (int j1 = 0; j1 < degree; j1++){
          unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
          total_payload = total_payload + opt->edges[edge_idx2].a;
          int node_idx2 = opt->edges[edge_idx2].t;
          int new_idx = opt->map2[node_idx2];
          //if (hyperedge_idx == 752586) printf("i1 = %d, i2 = %d, sep_idx = %d, newidx = %d, node_idx = %d, hyperedge_idx = %d \n", i1, i2, sep_idx, new_idx, node_idx2, hyperedge_idx);
          if (new_idx > sep_idx){ // edge contain vertex in V\B. 
            mark_payload = mark_payload + opt->edges[edge_idx2].a;
            //mark_edges[num_edges] = edge_idx;
            num_edges++;
          }
        }
        //if (hyperedge_idx == 752586) printf("normal mark, %d\n", num_edges);
        if (num_edges > 0 && num_edges < degree){
          //double weight = opt->hyperedges[i].weight / num_edges; //later change to weighted version. 
          flag_k = 1;
          ratio = total_payload / mark_payload;
          // printf("Ratio is %f\n", ratio);
          // printf("num edges %d, degree = %d\n", num_edges, degree);
          for (int j1 = 0; j1 < degree; j1++){
            unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
            int node_idx2 = opt->edges[edge_idx2].t;
            int new_idx = opt->map2[node_idx2];
            //nodes.push_back(opt->hypernodes[node_idx2]);
            if (new_idx > sep_idx){ // edge contain vertex in V\B. 
              //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->edges[edge_idx].a + weight;
              opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r + opt->edges[edge_idx2].a * ratio;
              //temp_a[edge_idx].a = opt->edges[edge_idx].a * ratio;
            }
          }
          //auto it = mySet.find(hyperedge_idx);
          if (mySet.count(hyperedge_idx) > 0) { //already inside mySet
            sep_idx = i1-1;
            int num_edges = 0;
            int degree = opt->hyperedges[hyperedge_idx].degree;
            //int *mark_edges = (int*)malloc(degree*sizeof(int));
            double ratio = 0.0;
            double total_payload = 0.0;
            double mark_payload = 0.0;
            for (int j1 = 0; j1 < degree; j1++){
              unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
              total_payload = total_payload + opt->edges[edge_idx2].a;
              int node_idx2 = opt->edges[edge_idx2].t;
              int new_idx = opt->map2[node_idx2];
              //if (new_idx >= i1 && new_idx <= sep_idx) printf("strange2\n");
              if (new_idx > sep_idx){ // edge contain vertex in V\B. 
                mark_payload = mark_payload + opt->edges[edge_idx2].a;
                //mark_edges[num_edges] = edge_idx;
                num_edges++;
              }
            }
            if (num_edges > 0 && num_edges < degree){
              //double weight = opt->hyperedges[i].weight / num_edges; //later change to weighted version. 
              flag_k = 1;
              ratio = total_payload / mark_payload;
              // printf("Ratio is %f\n", ratio);
              // printf("num edges %d, degree = %d\n", num_edges, degree);
              for (int j1 = 0; j1 < degree; j1++){
                unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
                int node_idx2 = opt->edges[edge_idx2].t;
                int new_idx = opt->map2[node_idx2];
                if (new_idx > sep_idx){ // edge contain vertex in V\B. 
                  //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->edges[edge_idx].a + weight;
                  opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r - opt->edges[edge_idx2].a * ratio;
                  //temp_a[edge_idx].a = opt->edges[edge_idx].a * ratio;
                }
                // else{
                //   //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->links[edge_idx].a + weight;
                //   opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r + opt->edges[edge_idx2].a;
                //   //temp_a[edge_idx].a = 0.0;
                // }
              }
            }
            sep_idx = i2;
          } else { // not found, then insert it into mySet
            if (hyperedge_idx == 1001865) printf("inserted\n");
            mySet.insert(hyperedge_idx);
            for (int j1 = 0; j1 < degree; j1++){
              unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
              int node_idx2 = opt->edges[edge_idx2].t;
              int new_idx = opt->map2[node_idx2];
              int flag_mark = 0;
              // if (new_idx < i1){
              //    printf("i1 = %d, i2 = %d, newidx = %d, node_idx = %d, hyperedge_idx = %d ", i1, i2, new_idx, node_idx2, hyperedge_idx);
              //    for (int j=0; j < opt->hypernodes[node_idx2].degree; j++){
              //       int edge_idx0 = opt->hypernodes[node_idx].neighboredges[j];
              //       int hyperedge_idx0 = opt->edges[edge_idx].s;
              //       if (hyperedge_idx0 == hyperedge_idx) flag_mark = 1;
              //    }
              //    if (flag_mark == 0) printf("strange\n");
              //    return 1;
              // }
              //nodes.push_back(opt->hypernodes[node_idx2]);
              //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->links[edge_idx].a + weight;
              opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r - opt->edges[edge_idx2].a;
              //temp_a[edge_idx].a = 0.0;
            }
            //printf("Next\n");
          }
        }
        else if (num_edges == 0){
          //mark_hyper[hyperedge_idx] = 1;
          //mySet.erase(hyperedge_idx);
          size_t numErased = mySet.erase(hyperedge_idx);
          
          //mark_hyper = 1;
          //auto it = mySet.find(hyperedge_idx);
          //numErased = 2;
          if (numErased>0) {
            sep_idx = i1-1;
            int num_edges = 0;
            int degree = opt->hyperedges[hyperedge_idx].degree;
            //int *mark_edges = (int*)malloc(degree*sizeof(int));
            double ratio = 0.0;
            double total_payload = 0.0;
            double mark_payload = 0.0;
            for (int j1 = 0; j1 < degree; j1++){
              unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
              total_payload = total_payload + opt->edges[edge_idx2].a;
              int node_idx2 = opt->edges[edge_idx2].t;
              int new_idx = opt->map2[node_idx2];
              //if (new_idx > i2) printf("strange2\n");
              if (new_idx > sep_idx){ // edge contain vertex in V\B. 
                mark_payload = mark_payload + opt->edges[edge_idx2].a;
                //mark_edges[num_edges] = edge_idx;
                num_edges++;
              }
            }
            if (num_edges > 0 && num_edges < degree){
              //double weight = opt->hyperedges[i].weight / num_edges; //later change to weighted version. 
              flag_k = 1;
              ratio = total_payload / mark_payload;
              // printf("Ratio is %f\n", ratio);
              // printf("num edges %d, degree = %d\n", num_edges, degree);
              for (int j1 = 0; j1 < degree; j1++){
                unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
                int node_idx2 = opt->edges[edge_idx2].t;
                int new_idx = opt->map2[node_idx2];
                if (new_idx > sep_idx){ // edge contain vertex in V\B. 
                  //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->edges[edge_idx].a + weight;
                  opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r - opt->edges[edge_idx2].a * ratio + opt->edges[edge_idx2].a;
                  //temp_a[edge_idx].a = opt->edges[edge_idx].a * ratio;
                }
                else{
                  //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->links[edge_idx].a + weight;
                  opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r + opt->edges[edge_idx2].a;
                  //temp_a[edge_idx].a = 0.0;
                }
              }
            }
          }
        }
      
      }
    }
    sep_idx = i2;
    double min_r = -1.0;
    double max_r = -1.0;
    for (unsigned j1 = 0; j1 < opt->n; j1++){
      unsigned node_idx = node_arr[j1].id;
      //printf("sep_idx = %d\n", sep_idx);
      if (j1 <= sep_idx){
        if (min_r < 0){
          min_r = opt->hypernodes[node_idx].r;
        }
        if (opt->hypernodes[node_idx].r< min_r){
          min_r = opt->hypernodes[node_idx].r;
        }
      }
      else{
        //printf("j1 is %d, hypernode rate is %f", j1, opt->hypernodes[node_idx].r);
        if (opt->hypernodes[node_idx].r > max_r){
          max_r = opt->hypernodes[node_idx].r;
          //printf("max r is %f", max_r);
        }
      }
    }
    if (min_r - max_r> 0){
      printf("min_r = %f, max_r = %f \n", min_r, max_r);
      decom[num++] = acc;
      printf("flagk = %d\n", flag_k);
      printf("Stable, %d\n", acc);
      //fprintf(file_arr, "%d, min_r = %f, max_r = %f \n",acc, min_r, max_r);
    }
    //fprintf(file_arr, "%d, min_r = %f, max_r = %f \n",acc, min_r, max_r);
    // for (const auto& n : nodes) {
    //   opt->hypernodes[n.id].r = temp_r[n.id];
    //   //std::cout << "Node ID: " << n.id << ", Weight: " << n.weight << std::endl;
    // }
    // nodes.clear();
    printf("One part %d checked\n", acc);
  }
  free(nag);
  free(mark_hyper);
  free(temp_r);
  //fclose(file_arr);
  //FILE *file_decomp = fopen("decomposition_mark.txt", "w");
  //for (int k = 0; k < num; k++) fprintf(file_decomp, "%d\n", decom[k]);
  //fclose(file_decomp);
  printf("Start to do verification using LP %d\n", num);

  // subgraph *sg = allocsubgraph(opt->n, opt->m, opt->e);
  
  // for (int i = 0; i < opt->n; i++){
  //   sg->hypernodes[i] = opt->hypernodes[i];
  //   sg->hypernodes[i].part = -1;
  // }
  // for (int k1 = 0; k1 < opt->m; k1++){
  //   sg->hyperedges[k1] = opt->hyperedges[k1];
  //   sg->hyperedges[k1].part = -1;
  // }
  // for (unsigned k1 = 0; k1 < opt->e; k1++){
  //   sg->edges[k1] = opt->edges[k1];
  //   sg->edges[k1].s = opt->edges[k1].s;
  //   sg->edges[k1].t = opt->edges[k1].t;
  // }
  for (int k = 0; k < num; k++){
    int i1 = 0;
    if (k > 0) i1 = decom[k-1];
    int i2 = decom[k]-1;
    printf("i1 = %d, i2 = %d \n", i1, i2);
    for (int k2=i1; k2 <=i2; k2++){
      int idx = node_arr[k2].id;
      opt->hypernodes[idx].part = (k+1);
    }
    // int n_nodes = i2 - i1+1;
    // printf("marb = %d, n_nodes = %d\n",marb - n_nodes, n_nodes);
    // marb = marb - n_nodes;
  }

  //freeopt(opt);
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  long memory_usage_mb = usage.ru_maxrss / 1024;
  printf("Memory usage: %ld MB\n", memory_usage_mb);
  //FILE *file_parts = fopen("parts_mark.txt", "w");
  int start_part = 0;
  unsigned removed_hyperedge = 0;
  unsigned removed_e = 0;
  // for (int i =0 ; i < opt->m; i++){
  //   int hyperedge_idx = i;
  //   for (int j2=0; j2 < opt->hyperedges[hyperedge_idx].degree; j2++){
  //     int link_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j2];
  //     printf("s = %d, hyperedge idx = %d\n", opt->edges[link_idx2].s, hyperedge_idx);
  //     int hypernode_idx0 = opt->edges[link_idx2].t;
  //     if (hyperedge_idx != opt->edges[link_idx2].s){
  //       exit(-1);
  //     }
  //   }
  // }
  int* included_mark = (int*)malloc(opt->m*sizeof(int));
  for (int i=0; i < opt->m;i++){
    included_mark[i]=0;
  }

  for (int k = start_part; k < num; k++){
    // if (k == num - 1){
    //   verify_maxflow(sg, Mark_nodes);
    //   free(decom);
    //   //freesubgraph(sg);
    //   fclose(file_parts);
    //   return 1;
    // }
    long memory_usage_kb = usage.ru_maxrss;
    printf("Start Memory usage: %ld KB\n", memory_usage_kb);
    int i1 = 0;
    if (k > 0) i1 = decom[k-1];
    int i2 = decom[k]-1;
    int n_nodes = i2 - i1+1;
    printf("Start to process part %d\n", k);
    subgraph *sg1=allocsubgraph(n_nodes,opt->m - removed_hyperedge, opt->e - removed_e);
    //subgraph *sg2=allocsubgraph(sg->n-n_nodes,sg->m, sg->e);
    //current how large
    sg1->n = n_nodes;
    printf("n_nodes = %d\n", n_nodes);
    unsigned j_edge = 0;
    unsigned j_hyperedge = 0;
    //map
    printf("Start to assign edges\n");
    int **temp_ptr = (int**)malloc((opt->m - removed_hyperedge)*sizeof(int*));
    for (int k2=i1; k2 <=i2; k2++){
      int idx = node_arr[k2].id;
      //printf("mark1\n");
      //printf("mark, k2 - i1 = %d, idx = %d\n", k2- i1, idx);
      sg1->hypernodes[k2-i1].id = k2-i1;
      sg1->hypernodes[k2-i1].old_id = opt->hypernodes[idx].old_id;
      sg1->hypernodes[k2-i1].weight = opt->hypernodes[idx].weight;
      //printf("k2 - i1 = %d, opt map is %d, idx is %d\n", k2-i1, opt->map2[idx], idx);
      //= opt->hypernodes[idx];
      //printf("degree %d\n", opt->hypernodes[idx].degree);
      //printf("mark, n_nodes = %d\n", n_nodes);
      //printf("%d\n", max_degree);

      for (int j =0; j < opt->hypernodes[idx].degree; j++){
        int link_idx = opt->hypernodes[idx].neighboredges[j];
        //printf("mark1.5, %d\n", link_idx);
        //printf("%d\n", max_degree);
        int hyperedge_idx = opt->edges[link_idx].s;
        //printf("mark2, %d\n", hyperedge_idx);
        int flag_include = 1;

        for (int j2=0; j2 < opt->hyperedges[hyperedge_idx].degree; j2++){
          int link_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j2];
          //printf("s = %d, mark3, %d\n", opt->edges[link_idx2].s, opt->hyperedges[hyperedge_idx].degree);
          //printf("j_hyperedge = %d, \n", j_hyperedge);
          int hypernode_idx0 = opt->edges[link_idx2].t;
          if (opt->hypernodes[hypernode_idx0].part > (k+1)) flag_include = -1;
          //if (opt->hypernodes[hypernode_idx0].part == (k+1) && flag_include>=0) flag_include = 1;
        }
        //printf("mark4\n");
        if (flag_include == 1 && included_mark[hyperedge_idx] == 0){ //应该包含在里面
          sg1->hyperedges[j_hyperedge] = opt->hyperedges[hyperedge_idx];
          included_mark[hyperedge_idx] = 1;
          sg1->hyperedges[j_hyperedge].neighboredges = (int*)malloc(opt->hyperedges[hyperedge_idx].degree * sizeof(int));
          temp_ptr[j_hyperedge] = sg1->hyperedges[j_hyperedge].neighboredges; // to free this later. 
          int j_temp = 0;
          for (int j2=0; j2 < opt->hyperedges[hyperedge_idx].degree; j2++){
            int link_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j2];
            int hypernode_idx0 = opt->edges[link_idx2].t;
            if (opt->hypernodes[hypernode_idx0].part == (k+1)){
              sg1->hyperedges[j_hyperedge].neighboredges[j_temp] = j_edge;
              sg1->edges[j_edge].s = j_hyperedge;
              sg1->edges[j_edge].t = opt->map2[hypernode_idx0] - i1;
              j_temp++;
              j_edge++;
            }
          }
          sg1->hyperedges[j_hyperedge].degree = j_temp;
          j_hyperedge++;
        }
      }
      //sg->hypernodes[idx].part = (k+1);
    }
    sg1->m = j_hyperedge;
    sg1->e = j_edge;
    printf("count = %d\n", count);
    removed_hyperedge = removed_hyperedge + sg1->m;
    removed_e = removed_e + sg1->e;
    printf("sg1 created\n");
    printf("sg1->m = %d, sg1->n = %d\n", sg1->m, sg1->n);
    //fprintf(file_parts, "sg1->m = %d, sg1->n = %d\n", sg1->m, sg1->n);
    //free(isin);
    //freesubgraph(sg);
    printf("Verify using maxflow\n");
    verify_maxflow(sg1, Mark_nodes);
    //printf("start to free sg2, opt- removed = %d, sg1->m= %d\n", opt->m - removed_hyperedge+ sg1->m, sg1->m);
    // printf("size of = %lu\n", sizeof(temp_ptr));
    // printf("size of = %lf\n", sizeof(temp_ptr) / sizeof(int*));
    for (int i =0; i < j_hyperedge; i++){
      free(temp_ptr[i]);
    }
    free(temp_ptr);
    //freesubgraph2(sg2);
    memory_usage_kb = usage.ru_maxrss;
    printf("End Memory usage: %ld KB\n", memory_usage_kb);
  }
  //fclose(file_parts);
  return 1;
}






int main(int argc,char** argv){
	optim* opt;
  //isoreg *fit,*fit2;
	unsigned nsgs,i,k;
	unsigned nthreads=4;
	unsigned rep=150;
	double *err;
  char* edgelist="bipartite_mark.txt";
  // char* rates=argv[4];
  // char* pavafit=argv[5];
	// char* cuts=argv[6];

  //omp_set_num_threads(nthreads);

	time_t t0,t1,t2,t3, t4;
	t1=time(NULL);
	t0=t1;
  printf("- Reading edgelist from file %s\n",edgelist);
	opt=readedgelist(edgelist);
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  long memory_usage_mb = usage.ru_maxrss / 1024;
  printf("Memory usage: %ld MB\n", memory_usage_mb);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Building the datastructure\n");
	t1=time(NULL);
	relabel(opt); 
  getrusage(RUSAGE_SELF, &usage);
  memory_usage_mb = usage.ru_maxrss / 1024 / 1024;
  printf("Memory usage: %ld GB\n", memory_usage_mb);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Building the datastructure\n");
	printf("- Number of nodes = %u\n",opt->n);
	printf("- Number of edges = %llu\n",opt->e);
	printf("- Computing the locally densest decomposition\n");
  printf("- Step 1: Frank-Wolf gradiant descent (%u iterations)\n",rep);
	t1=time(NULL);
	init(opt);
  getrusage(RUSAGE_SELF, &usage);
  memory_usage_mb = usage.ru_maxrss / 1024 / 1024;
  printf("Memory usage: %ld GB\n", memory_usage_mb);
  //read weights
  //read_weights(opt, "bipartite_mark_E.txt", "bipartite_mark_N.txt");
  // if we use FISTA, we need the following. Otherwise, we can comment them. 
  double max_degree = get_max_weighted_degree(opt);
	link0 *z; 
	z = (link0*)malloc(opt->e*sizeof(link0));
	link0 *y; 
	y = (link0*)malloc(opt->e*sizeof(link0));
	link0 *last_x;
	last_x = (link0*)malloc(opt->e*sizeof(link0));
  getrusage(RUSAGE_SELF, &usage);
  memory_usage_mb = usage.ru_maxrss / 1024 / 1024;
  printf("Memory usage: %ld GB\n", memory_usage_mb);
	for (i = 0; i < opt->e; i++){
		z[i].s = opt->edges[i].s;
		z[i].t = opt->edges[i].t;
		z[i].a = opt->edges[i].a;
		//z[i].b = 1 - opt->edges[i].a;

		y[i].s = opt->edges[i].s;
		y[i].t = opt->edges[i].t;
		y[i].a = opt->edges[i].a;

		last_x[i].s = opt->edges[i].s;
		last_x[i].t = opt->edges[i].t;
		last_x[i].a = opt->edges[i].a;
		//last_x[i].b = 1 - opt->edges[i].a;
	}
  printf("start the algorithm...\n");
  //Out_put_density(opt, "output_init.txt");
  //---------------------------------------
  double lambda = 0;
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, &lambda);
		//onepass_PR(opt);
    onepass_PR_plus_product_normal_2(opt, last_x, z);
    printf("iteration %d finished\n", i+1);
		//onepass_FW1(opt);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
	}
  
  update_density(opt);
	t2=time(NULL);
    printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
  //Out_put_density(opt, "output_PR1.txt");
  //Out_put_weights(opt, "output_PR1w.txt");
  printf("%d, %d\n", opt->m, opt->n);
  node *Mark_nodes = (node*)malloc(opt->n * sizeof(node));
  for (int i =0; i < opt->n; i++){
    Mark_nodes[i].old_id = opt->hypernodes[i].old_id;
    Mark_nodes[i].r = opt->hypernodes[i].r;
    Mark_nodes[i].weight = opt->hypernodes[i].weight;
  }
  //qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes);

  t3=time(NULL);
  int res = verify(opt, Mark_nodes);
  printf("verify, %d", res);
  t4 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n",(t4-t3)/3600,((t4-t3)%3600)/60,((t4-t3)%60));

  FILE *exact = fopen("Exact_normal.txt", "w");

  for (int i = 0; i < opt->n; i++){
    fprintf(exact, "%d %.16f\n", i, Mark_nodes[i].r_correct);
    if (Mark_nodes[i].old_id != i) printf("strange");
  }
  fclose(exact);

}



/*
    unsigned *newlabel=(unsigned*)malloc(sg->n*sizeof(unsigned));
    unsigned j1=0;
    for(unsigned i=0;i<sg->n;i++){
      if (sg->hypernodes[i].part ==(k+1)){
        sg1->hypernodes[j1] = sg->hypernodes[i];
        newlabel[i]=j1++;
        //printf("One part i = %d, j1 = %d, n_nodes = %d\n", i, j1, n_nodes);
      }
      else{
        //printf("Part = %d\n", sg->hypernodes[i].part);
        sg2->hypernodes[j2] = sg->hypernodes[i];
        newlabel[i]=j2++;
      }
    }
    sg1->n = j1;
    sg2->n = j2;
    unsigned j1_edge = 0;
    unsigned j2_edge = 0;
    j1=0, j2=0;
    printf("Start to assign edges\n");
    for (int k2=i1; k2 <=i2; k2++){
      for ()
      sg->hypernodes[idx].part = (k+1);
    }

    for (int i = 0; i < sg->m; i++){
      int flag = 0;
      for (int j = 0; j < sg->hyperedges[i].degree; j++){
        int edge_idx = sg->hyperedges[i].neighboredges[j];
        if (sg->hypernodes[sg->edges[edge_idx].t].part !=(k+1)){
          flag = 1;
        }
      }
      if (flag == 0){
        sg1->hyperedges[j1] = sg->hyperedges[i];
        for (int j = 0; j < sg->hyperedges[i].degree; j++){
          int edge_idx = sg->hyperedges[i].neighboredges[j];
          sg1->edges[j1_edge].s = j1;
          sg1->edges[j1_edge].t = newlabel[sg->edges[edge_idx].t];
          sg1->hyperedges[j1].neighboredges[j] = j1_edge;
          j1_edge++;
        }
        j1++;
      }
      else{
        sg2->hyperedges[j2] = sg->hyperedges[i];
        int j_tmp = 0;
        for (int j = 0; j < sg->hyperedges[i].degree; j++){
          int edge_idx = sg->hyperedges[i].neighboredges[j];
          if (sg->hypernodes[sg->edges[edge_idx].t].part !=(k+1)){
            sg2->edges[j2_edge].s = j2;
            sg2->edges[j2_edge].t = newlabel[sg->edges[edge_idx].t];
            sg2->hyperedges[j2].neighboredges[j_tmp] = j2_edge;
            j2_edge++;
            j_tmp++;
          }
        }
        sg2->hyperedges[j2].degree = j_tmp; 
        j2++;
      }
    }
    sg1->m = j1;
    sg2->m = j2;
    sg1->e=j1_edge;
    sg2->e=j2_edge;



    for (const int& hyperedge_idx : mySet) {
      if (mark_hyper[hyperedge_idx] == 0){
        int num_edges = 0;
        int degree = opt->hyperedges[hyperedge_idx].degree;
        //int *mark_edges = (int*)malloc(degree*sizeof(int));
        double ratio = 0.0;
        double total_payload = 0.0;
        double mark_payload = 0.0;
        for (int j1 = 0; j1 < degree; j1++){
          unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
          total_payload = total_payload + opt->edges[edge_idx2].a;
          int node_idx2 = opt->edges[edge_idx2].t;
          int new_idx = opt->map2[node_idx2];
          if (new_idx > sep_idx){ // edge contain vertex in V\B. 
            mark_payload = mark_payload + opt->edges[edge_idx2].a;
            //mark_edges[num_edges] = edge_idx;
            num_edges++;
          }
        }
        if (num_edges > 0 && num_edges < degree){
          //double weight = opt->hyperedges[i].weight / num_edges; //later change to weighted version. 
          flag_k = 1;
          ratio = total_payload / mark_payload;
          // printf("Ratio is %f\n", ratio);
          // printf("num edges %d, degree = %d\n", num_edges, degree);
          for (int j1 = 0; j1 < degree; j1++){
            unsigned edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j1];
            int node_idx2 = opt->edges[edge_idx2].t;
            int new_idx = opt->map2[node_idx2];
            nodes.push_back(opt->hypernodes[node_idx2]);
            if (new_idx > sep_idx){ // edge contain vertex in V\B. 
              //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->edges[edge_idx].a + weight;
              opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r + opt->edges[edge_idx2].a * ratio - opt->edges[edge_idx2].a;
              //temp_a[edge_idx].a = opt->edges[edge_idx].a * ratio;
            }
            else{
              //opt->hypernodes[node_idx].r =opt->hypernodes[node_idx].r - opt->links[edge_idx].a + weight;
              opt->hypernodes[node_idx2].r = opt->hypernodes[node_idx2].r - opt->edges[edge_idx2].a;
              //temp_a[edge_idx].a = 0.0;
            }
          }
        }
        if (num_edges == 0){
          mark_hyper[hyperedge_idx] = 1;
          mySet.erase(hyperedge_idx);
          // auto it = mySet.find(hyperedge_idx);
          // if (it != mySet.end()) {
          // }
        }
      }
    }
*/