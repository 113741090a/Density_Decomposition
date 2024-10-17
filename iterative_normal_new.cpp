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
#include <sys/time.h>
#include <math.h>

// #include <stdio.h>
// #include <stdlib.h>
#include <omp.h>



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
  //used for greedy++
  int heap_idx;
  double d_G;
  //-------------------
    double weight;
	  double r;//rate value of node id
    double r_correct; //correct value of the rate value;
    double rho; //density of node id. r/weight.
    int *neighboredges;
    int degree;
    int part;
    //node *neighbors;
    // 重载小于运算符，用于优先队列比较
    bool operator<(const node& other) const {
        return r < other.r;
    }
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
  double *ne;//ne[i]=number of edges from i to nodes before (used for pava)
	unsigned *cd;//cumulative degree
	unsigned *cuts;
  unsigned iter;//number of iterations
  double tk;
} optim;

void Out_put_density(optim *opt, char *file_name);
void Out_put_weights(optim *opt, char *file_name, int rep);


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
			opt->edges=(link0*)realloc(opt->edges,e1*sizeof(link0) );
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
	for (k=0;k<opt->e;k++){
    int hyperedge_idx = opt->edges[k].s;
    double weight = opt->hyperedges[hyperedge_idx].weight;
    int deg = opt->hyperedges[hyperedge_idx].degree;
    //opt->edges[k].a= weight / 2.0;
    opt->edges[k].a= weight / ((double) deg);

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
    int deg = opt->hyperedges[hyperedge_idx].degree;
    //opt->edges[k].a= weight / 2.0; //already updated
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
  //calculate d_G
  for (k=0;k<opt->n;k++){
    for (int x=0; x < opt->hypernodes[k].degree; x++){
      int edge_idx = opt->hypernodes[k].neighboredges[x];
      int hyperedge_idx = opt->edges[edge_idx].s;
      double weight = opt->hyperedges[hyperedge_idx].weight;
      opt->hypernodes[k].d_G = opt->hypernodes[k].d_G + weight;
    }
  }
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


void quicksort(int arr[], int low, int high, link0* edges);

void update_density2(optim *opt, link0* weight);

void update_density(optim *opt);

void onepass_PR_plus_product_normal(optim *opt, link0* last_x, link0* z){
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
    if (isnan(z[k].a)){
      printf("sgrange\n");
      exit(1);
    }
    //exp((1+gamma)*log(opt->edges[k].a) - gamma*log(last_x[k].a));
  }


  //do projection. In a more easier way. 
  //#pragma omp parallel for private(i,j,k)
  for (i=0; i < opt->m; i++){
    //first step, sort neighboring edges of i, according to y_{ij}
    //qsort(opt->hyperedges[i].neighboredges,opt->hyperedges[i].degree,sizeof(edge),compare_edges);
    //printf("\n");
    int edge_idx1 = opt->hyperedges[i].neighboredges[0];
    int edge_idx2 = opt->hyperedges[i].neighboredges[1];
    double hypernode_weight = opt->hyperedges[i].weight;
    double edge_weight1 = z[edge_idx1].a;
    double edge_weight2 = z[edge_idx2].a;
    z[edge_idx1].a = (hypernode_weight + edge_weight1 - edge_weight2)/(2.0);
    if (z[edge_idx1].a < 0) z[edge_idx1].a = 0.0;
    if (z[edge_idx1].a > hypernode_weight) z[edge_idx1].a = hypernode_weight;
    z[edge_idx2].a = hypernode_weight - z[edge_idx1].a;
    //printf("\n");
  }
  // for (k = 0; k < opt->e; k++){
  //   if (z[k].a < 0 || z[k].a > 1) printf("strange\n");
  //   //exp((1+gamma)*log(opt->edges[k].a) - gamma*log(last_x[k].a));
  // }

  // update_density(opt);
  // double res1 = 0.0;
  // for (int i = 0; i < opt->n; i++){
  //   res1 = res1 + opt->hypernodes[i].r * opt->hypernodes[i].r;
  // }
  // update_density2(opt, last_x);
  // double res2 = 0.0;
  // for (int i = 0; i < opt->n; i++){
  //   res2 = res2 + opt->hypernodes[i].r * opt->hypernodes[i].r;
  // }
  // printf("res1 = %lf, res2 = %lf\n", res1, res2);
  //*lambda = lambda1;
}



//proportional projection
void onepass_PR_plus_product_normal_pro(optim *opt, link0* last_x, link0* z){
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
    //first step, sort neighboring edges of i, according to y_{ij}
    //qsort(opt->hyperedges[i].neighboredges,opt->hyperedges[i].degree,sizeof(edge),compare_edges);
    //printf("\n");

    int edge_idx1 = opt->hyperedges[i].neighboredges[0];
    int edge_idx2 = opt->hyperedges[i].neighboredges[1];
    double hypernode_weight = opt->hyperedges[i].weight;
    double edge_weight1 = z[edge_idx1].a;
    double edge_weight2 = z[edge_idx2].a;
    double parameter = hypernode_weight / (edge_weight1 + edge_weight2);
    if (parameter != 1.0){
      z[edge_idx1].a = parameter * z[edge_idx1].a;
      z[edge_idx2].a = parameter * z[edge_idx2].a;
    }
    //printf("\n");
  }
}


void onepass_PR_plus_lin_normal(optim *opt, link0* last_x, link0* z){
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
    int mark_flag = 0;
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double r = opt->hypernodes[neighbor_node_idx].r;
      if (r == 0.0){
        mark_flag = 1;
        payload = payload + opt->hypernodes[neighbor_node_idx].weight * 1 / opt->hypernodes[neighbor_node_idx].degree; 
      }
      else payload = payload + opt->hypernodes[neighbor_node_idx].weight * z[edge_idx].a / r;
    }
//be carefull here. 
    if (mark_flag == 1){
      int count = 0;
      for (j=0; j < opt->hyperedges[i].degree; j++){
        int edge_idx = opt->hyperedges[i].neighboredges[j];
        int neighbor_node_idx = opt->edges[edge_idx].t;
        double r = opt->hypernodes[neighbor_node_idx].r;
        if (r == 0.0){
          opt->edges[edge_idx].a = opt->hyperedges[i].weight * (opt->hypernodes[neighbor_node_idx].weight * 1 / opt->hypernodes[neighbor_node_idx].degree)/payload;
        }
        else{
          opt->edges[edge_idx].a = opt->hyperedges[i].weight * (opt->hypernodes[neighbor_node_idx].weight * z[edge_idx].a / r)/payload;
        }
      }
      continue;
    }
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
      if (isnan(opt->edges[edge_idx].a)){
        printf("sgrange, r = %lf, payload = %lf,\n", r, payload);
        exit(1);
      }
    }
  }

  //linear combination. 
  //double lambda1 = (1 + sqrt(1 + 4 * pow(*lambda, 2))) /2.0;
  //double gamma = (double)(t-1)/(double)(t+2);
  double gamma = (double)(t-1)/(double)(t+2);
  for (k = 0; k < opt->e; k++){
    z[k].a = opt->edges[k].a + gamma * (opt->edges[k].a - last_x[k].a);
    //printf("values %lf, %lf\n", opt->edges[k].a, last_x[k].a);
    if (isnan(z[k].a)){
      printf("sgrange\n");
      exit(1);
    }
  }



  //do projection. In a more easier way. 
  //#pragma omp parallel for private(i,j,k)
  for (i=0; i < opt->m; i++){
    //first step, sort neighboring edges of i, according to y_{ij}
    //qsort(opt->hyperedges[i].neighboredges,opt->hyperedges[i].degree,sizeof(edge),compare_edges);
    //printf("\n");
    int edge_idx1 = opt->hyperedges[i].neighboredges[0];
    int edge_idx2 = opt->hyperedges[i].neighboredges[1];
    double edge_weight1 = z[edge_idx1].a;
    double edge_weight2 = z[edge_idx2].a;
    double hypernode_weight = opt->hyperedges[i].weight;
    z[edge_idx1].a = (hypernode_weight + edge_weight1 - edge_weight2)/(2.0);
    if (z[edge_idx1].a < 0) z[edge_idx1].a = 0.0;
    if (z[edge_idx1].a > hypernode_weight) z[edge_idx1].a = hypernode_weight;
    z[edge_idx2].a = hypernode_weight - z[edge_idx1].a;
    //printf("\n");
  }
  //*lambda = lambda1;
}



void update_density(optim *opt);






//one pass over all edges, FW, sequential, don't use it. 
void onepass_FW(optim *opt){
	unsigned i,j;
	unsigned long long k;
	double gamma;
	opt->iter++;
	gamma=2./(2.+opt->iter);
  for (i=0;i<opt->m;i++){
    int min_idx = 0;
    //--in order to initialize min_r-----:
    int edge_idx = opt->hyperedges[i].neighboredges[0];
    int neighbor_node_idx = opt->edges[edge_idx].t;
    double min_r = opt->hypernodes[neighbor_node_idx].r / opt->hypernodes[neighbor_node_idx].weight;

    for (j=1; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double r = opt->hypernodes[neighbor_node_idx].r  / opt->hypernodes[neighbor_node_idx].weight;;
      if (r < min_r){
        min_r = r;
        min_idx = j;
      }
    }
//be carefull here. 
    double hyperedge_weight = opt->hyperedges[i].weight;
    for (j=0; j < opt->hyperedges[i].degree; j++){
      if (j == min_idx){
        int edge_idx = opt->hyperedges[i].neighboredges[j];
        int neighbor_node_idx = opt->edges[edge_idx].t;
        //opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * (1 - opt->edges[edge_idx].a); //weight
        opt->edges[edge_idx].a = (1. - gamma)*opt->edges[edge_idx].a + gamma * hyperedge_weight; // 1 is the weight of hyperedge. 
      }
      else{
        int edge_idx = opt->hyperedges[i].neighboredges[j];
        int neighbor_node_idx = opt->edges[edge_idx].t;
        //opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r - gamma * (opt->edges[edge_idx].a); //weight
        opt->edges[edge_idx].a = (1. - gamma)*opt->edges[edge_idx].a; // 1 is the weight of hyperedge. 
      }
    }
  }
  update_density(opt);
}


//one pass over all edges
void onepass_FW_syn(optim *opt){
	unsigned i,j;
	unsigned long long k;
	opt->iter++;
  double gamma = 2.0/(2.0+opt->iter);
  double payload = 0;
  for (i =0; i < opt->e; i++){
    opt->edges[i].a = (1 - gamma) * opt->edges[i].a;
  }
  update_density(opt);
	//#pragma omp parallel for private(i,j,k)
  for (i=0;i<opt->m;i++){
    int min_idx = 0;
    int edge_idx = opt->hyperedges[i].neighboredges[min_idx];
    int neighbor_node_idx = opt->edges[edge_idx].t;
    double min_rho = opt->hypernodes[neighbor_node_idx].r / opt->hypernodes[neighbor_node_idx].weight;
    for (j=1; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double rho = opt->hypernodes[neighbor_node_idx].r / opt->hypernodes[neighbor_node_idx].weight;
      if (rho < min_rho){
        min_idx = j;
        min_rho = rho;
      }
    }
//be carefull here. 
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double weight = opt->hyperedges[i].weight;
      if (j == min_idx){
        //opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * (weight); 
        opt->edges[edge_idx].a = opt->edges[edge_idx].a + gamma * weight; //update r here 
        
      }
      // else{
      //   opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * (-opt->edges[edge_idx].a); //update r here, maybe need to change
      //   opt->edges[edge_idx].a = (1.0 - gamma)* opt->edges[edge_idx].a;
      // }
    }
  }
}



//quicksort, neighboring edges of i, according to the y_{ij}
void swap2(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int partition2(int arr[], int low, int high, optim *opt) {
    int pivot_idx = arr[high];
    int hypernode_idx = opt->edges[pivot_idx].t;
    double pivot_value = opt->hypernodes[hypernode_idx].r / opt->hypernodes[hypernode_idx].weight;
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
      int hypernode_idx2 = opt->edges[arr[j]].t;
      double j_value = opt->hypernodes[hypernode_idx2].r / opt->hypernodes[hypernode_idx2].weight;
      if (j_value < pivot_value) {
          i++;
          swap2(&arr[i], &arr[j]);
      }
    }
    swap2(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quicksort2(int arr[], int low, int high, optim *opt) {
    if (low < high) {
        int pi = partition2(arr, low, high, opt);
        quicksort2(arr, low, pi - 1, opt);
        quicksort2(arr, pi + 1, high, opt);
    }
}


//one pass over all edges
void onepass_FW_Water(optim *opt){
	unsigned i,j;
	unsigned long long k;
	opt->iter++;
  double gamma = 2.0/(2.0+opt->iter);
  for (i =0; i < opt->e; i++){
    opt->edges[i].a = (1 - gamma) * opt->edges[i].a;
  }
  update_density(opt);

	//#pragma omp parallel for private(i,j,k)
  for (i=0;i<opt->m;i++){
    quicksort2(opt->hyperedges[i].neighboredges, 0, opt->hyperedges[i].degree-1, opt);
    double mark = 0.0;
    double total_weight = 0.0;
    //printf("i = %d\n", i);
    // int min_idx = 0;
    // int edge_idx = opt->hyperedges[i].neighboredges[min_idx];
    // int neighbor_node_idx = opt->edges[edge_idx].t;
    // double min_rho = opt->hypernodes[neighbor_node_idx].r / opt->hypernodes[neighbor_node_idx].weight;
    double hyperedge_weight = opt->hyperedges[i].weight;
    //printf("degree = %d\n", opt->hyperedges[i].degree);
    if (opt->hyperedges[i].degree != 2){
      printf("strange things happens\n");
    }
    if (opt->hyperedges[i].degree == 1){
      int edge_idx = opt->hyperedges[i].neighboredges[0];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      opt->edges[edge_idx].a = opt->edges[edge_idx].a + gamma * hyperedge_weight;
      opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * hyperedge_weight;
      continue;
    }
    for (j=1; j < opt->hyperedges[i].degree; j++){
      //printf("enter the loop\n");
      int edge_idx1 = opt->hyperedges[i].neighboredges[j-1];
      int neighbor_node_idx1 = opt->edges[edge_idx1].t;
      double rho1 = opt->hypernodes[neighbor_node_idx1].r / opt->hypernodes[neighbor_node_idx1].weight;
      int edge_idx2 = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx2 = opt->edges[edge_idx2].t;
      double rho2 = opt->hypernodes[neighbor_node_idx2].r / opt->hypernodes[neighbor_node_idx2].weight;
      total_weight = total_weight + opt->hypernodes[neighbor_node_idx1].weight;
      //printf("%lf %lf, ", rho1, rho2);
      //printf("hyper edge weight is %lf\n", hyperedge_weight);
      if (mark + total_weight * (rho2 - rho1) / gamma < hyperedge_weight){
        mark = mark + total_weight * (rho2 - rho1) / gamma;
        //printf("mark is %lf\n", mark);
        double payload = 0.0;
        if (j == opt->hyperedges[i].degree - 1){
          total_weight = total_weight + opt->hypernodes[neighbor_node_idx2].weight;
          double rho = rho2 + gamma* (hyperedge_weight - mark) / total_weight;
          for (int j1=0; j1 <= j; j1++){
            //printf("makr\n");
            int edge_idx = opt->hyperedges[i].neighboredges[j1];
            int neighbor_node_idx = opt->edges[edge_idx].t;
            payload = payload + (rho * opt->hypernodes[neighbor_node_idx].weight - opt->hypernodes[neighbor_node_idx].r) / gamma;
            //printf("i = %d, idx = %d, payload = %lf\n",i, neighbor_node_idx, payload);
            opt->edges[edge_idx].a = opt->edges[edge_idx].a + (rho * opt->hypernodes[neighbor_node_idx].weight - opt->hypernodes[neighbor_node_idx].r);
            opt->hypernodes[neighbor_node_idx].r = rho * opt->hypernodes[neighbor_node_idx].weight;
          }
          //printf("good\n");
          if (abs(opt->hyperedges[i].weight- payload) > 0.0000001) printf("makr = %lf, weight = %lf, payload = %lf, rho1 = %lf, rho2 = %lf\n",mark, opt->hyperedges[i].weight, payload, rho1, rho2);
          payload = 0.0;
          for (int j1=0; j1 <= j; j1++){
            int edge_idx = opt->hyperedges[i].neighboredges[j1];
            payload = payload + opt->edges[edge_idx].a;
          }
          //printf("mark weight = %lf, payload = %lf\n", opt->hyperedges[i].weight, payload);
          if (abs(opt->hyperedges[i].weight- payload) > 0.0000001) exit(-1);
          break;
        }
      }
      else{
        double rho = rho1 + gamma*(hyperedge_weight - mark) / total_weight;
        double payload = 0.0;
        for (int j1=0; j1 < j; j1++){
          int edge_idx = opt->hyperedges[i].neighboredges[j1];
          int neighbor_node_idx = opt->edges[edge_idx].t;
          payload = payload + (rho * opt->hypernodes[neighbor_node_idx].weight - opt->hypernodes[neighbor_node_idx].r) / gamma;
          opt->edges[edge_idx].a = opt->edges[edge_idx].a + (rho * opt->hypernodes[neighbor_node_idx].weight - opt->hypernodes[neighbor_node_idx].r);
          opt->hypernodes[neighbor_node_idx].r = rho * opt->hypernodes[neighbor_node_idx].weight;
        }
        if (abs(opt->hyperedges[i].weight- payload) > 0.0000001) printf("123123weight = %lf, payload = %lf\n", opt->hyperedges[i].weight, payload);
        payload = 0.0;
        for (int j1=0; j1 < opt->hyperedges[i].degree; j1++){
          int edge_idx = opt->hyperedges[i].neighboredges[j1];
          payload = payload + opt->edges[edge_idx].a;
        }
        //printf("mark 123 weight = %lf, payload = %lf\n", opt->hyperedges[i].weight, payload);
        if (abs(opt->hyperedges[i].weight- payload) > 0.0000001) exit(-1);
        break;
      }
    }
    if (i == 1285){
      //printf("mark\n");
      double payload = 0.0;
      for (int j1=0; j1 < opt->hyperedges[i].degree; j1++){
        int edge_idx = opt->hyperedges[i].neighboredges[j1];
        payload = payload + opt->edges[edge_idx].a;
      }
      //printf("mark 123 weight = %lf, payload = %lf\n", opt->hyperedges[i].weight, payload);
      if (abs(opt->hyperedges[i].weight- payload) > 0.0000001) exit(-1);
    }
//be carefull here. 
      // else{
      //   opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * (-opt->edges[edge_idx].a); //update r here, maybe need to change
      //   opt->edges[edge_idx].a = (1.0 - gamma)* opt->edges[edge_idx].a;
      // }
  }
  for (i=0;i<opt->m;i++){
    double payload=0.0;
    for (j=0; j < opt->hyperedges[i].degree;j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      payload = payload + opt->edges[edge_idx].a;
    }
    if (abs(opt->hyperedges[i].weight- payload) > 0.0000001){
       printf("%d, weight = %lf, payload = %lf\n",i, opt->hyperedges[i].weight, payload);
       exit(-1);
    }
  }
}





//one pass over all edges
void onepass_FW_Elist(optim *opt){
	unsigned i,j;
	unsigned long long k;
	opt->iter++;
  double gamma = 2.0/(2.0+opt->iter);
  double payload = 0;
  for (i =0; i < opt->e; i++){
    opt->edges[i].a = (1 - gamma) * opt->edges[i].a;
  }
  update_density(opt);
	//#pragma omp parallel for private(i,j,k)
  for (i=0;i<opt->m;i++){
    int min_idx = 0;
    int edge_idx = opt->hyperedges[i].neighboredges[min_idx];
    int neighbor_node_idx = opt->edges[edge_idx].t;
    double min_rho = opt->hypernodes[neighbor_node_idx].r / opt->hypernodes[neighbor_node_idx].weight;
    for (j=1; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double rho = opt->hypernodes[neighbor_node_idx].r / opt->hypernodes[neighbor_node_idx].weight;
      if (rho < min_rho){
        min_idx = j;
        min_rho = rho;
      }
    }
//be carefull here. 
    for (j=0; j < opt->hyperedges[i].degree; j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int neighbor_node_idx = opt->edges[edge_idx].t;
      double weight = opt->hyperedges[i].weight;
      if (j == min_idx){
        opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * (weight); 
        opt->edges[edge_idx].a = opt->edges[edge_idx].a + gamma * weight; //update r here 
        
      }
      // else{
      //   opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + gamma * (-opt->edges[edge_idx].a); //update r here, maybe need to change
      //   opt->edges[edge_idx].a = (1.0 - gamma)* opt->edges[edge_idx].a;
      // }
    }
  }
  // for (i=0; i < opt->n; i++){
  //   opt->hypernodes[i].r = 0;
  // }
  // for (i=0; i < opt->e; i++){
  //   int hyper_node_idx = opt->edges[i].t;
  //   opt->hypernodes[hyper_node_idx].r = opt->hypernodes[hyper_node_idx].r + opt->edges[i].a;
  // }
  // for (i=0; i < opt->n; i++){
  //   opt->hypernodes[i].r = 0;
  // }
  // for (i=0;i<opt->m;i++){
  //   for (j=0; j < opt->hyperedges[i].degree; j++){
  //     int edge_idx = opt->hyperedges[i].neighboredges[j];
  //     int neighbor_node_idx = opt->edges[edge_idx].t;
  //     opt->hypernodes[neighbor_node_idx].r = opt->hypernodes[neighbor_node_idx].r + opt->edges[edge_idx].a;
  //   }
  // }
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

//quicksort, neighboring edges of i, according to the y_{ij}
void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(int arr[], int low, int high, link0* edges) {
    int pivot_idx = arr[high];
    double pivot_value = edges[pivot_idx].a;
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
      double j_value = edges[arr[j]].a;
      if (j_value < pivot_value) {
          i++;
          swap(&arr[i], &arr[j]);
      }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quicksort(int arr[], int low, int high, link0* edges) {
    if (low < high) {
        int pi = partition(arr, low, high, edges);

        quicksort(arr, low, pi - 1, edges);
        quicksort(arr, pi + 1, high, edges);
    }
}



//one pass over all edges
void onepass_FISTA_normal(optim *opt, double max_degree, link0* z, link0* y, link0* last_x){
	int i,j;
	unsigned long long k;
	opt->iter++;
  //in their own implementation, perhaps it is 0.9/Delta
  double alpha = 0.5 / (max_degree); //it's important to set 0.5 here. 
  double payload = 0;
  update_density2(opt, y);
  //#pragma omp parallel for private(i,j,k)
	for (k=0;k<opt->e;k++){//parfor
    int node_idx = z[k].t;
    z[k].a = y[k].a - 2.0*alpha*opt->hypernodes[node_idx].r / opt->hypernodes[node_idx].weight;
	}
  //do projection. In a more easier way. 
  //#pragma omp parallel for private(i,j,k)
  for (i=0; i < opt->m; i++){
    //first step, sort neighboring edges of i, according to y_{ij}
    //qsort(opt->hyperedges[i].neighboredges,opt->hyperedges[i].degree,sizeof(edge),compare_edges);
    //printf("\n");
    int edge_idx1 = opt->hyperedges[i].neighboredges[0];
    int edge_idx2 = opt->hyperedges[i].neighboredges[1];
    double hyperedge_weight = opt->hyperedges[i].weight;
    double edge_weight1 = z[edge_idx1].a;
    double edge_weight2 = z[edge_idx2].a;
    opt->edges[edge_idx1].a = (hyperedge_weight + edge_weight1 - edge_weight2)/(2.0);
    if (opt->edges[edge_idx1].a < 0) opt->edges[edge_idx1].a = 0;
    if (opt->edges[edge_idx1].a > hyperedge_weight) opt->edges[edge_idx1].a = hyperedge_weight;
    opt->edges[edge_idx2].a = hyperedge_weight - opt->edges[edge_idx1].a;
    //printf("\n");
  }
  for (k = 0; k < opt->e; k++){
    y[k].a = opt->edges[k].a + (opt->iter - 1.0) / (opt->iter + 2.0) * (opt->edges[k].a - last_x[k].a);
  }
  for (k = 0; k < opt->e; k++){
    last_x[k].a = opt->edges[k].a;
  }
}

//one pass over all edges
void onepass_FISTA_star_normal(optim *opt, double max_degree, link0* z, link0* y, link0* last_x){
	int i,j;
	unsigned long long k;
	opt->iter++;
  //in their own implementation, perhaps it is 0.9/Delta
  double alpha = 0.9 / (max_degree); //it's important to set 0.5 here. 
  double payload = 0;
  update_density2(opt, y);
  //#pragma omp parallel for private(i,j,k)
	for (k=0;k<opt->e;k++){//parfor
    int node_idx = z[k].t;
    z[k].a = y[k].a - 2.0*alpha*opt->hypernodes[node_idx].r / opt->hypernodes[node_idx].weight;
	}
  //do projection. In a more easier way. 
  //#pragma omp parallel for private(i,j,k)
  for (i=0; i < opt->m; i++){
    //first step, sort neighboring edges of i, according to y_{ij}
    //qsort(opt->hyperedges[i].neighboredges,opt->hyperedges[i].degree,sizeof(edge),compare_edges);
    //printf("\n");
    int edge_idx1 = opt->hyperedges[i].neighboredges[0];
    int edge_idx2 = opt->hyperedges[i].neighboredges[1];
    double hyperedge_weight = opt->hyperedges[i].weight;
    double edge_weight1 = z[edge_idx1].a;
    double edge_weight2 = z[edge_idx2].a;
    opt->edges[edge_idx1].a = (hyperedge_weight + edge_weight1 - edge_weight2)/(2.0);
    if (opt->edges[edge_idx1].a < 0) opt->edges[edge_idx1].a = 0;
    if (opt->edges[edge_idx1].a > hyperedge_weight) opt->edges[edge_idx1].a = hyperedge_weight;
    opt->edges[edge_idx2].a = hyperedge_weight - opt->edges[edge_idx1].a;
    //printf("\n");
  }
  for (k = 0; k < opt->e; k++){
    y[k].a = opt->edges[k].a + (opt->iter - 1.0) / (opt->iter + 2.0) * (opt->edges[k].a - last_x[k].a);
  }
  for (k = 0; k < opt->e; k++){
    last_x[k].a = opt->edges[k].a;
  }
}



void heap_swap(node* heap, int i1, int i2, optim* opt){
  node tmp = heap[i1];
  heap[i1] = heap[i2];
  heap[i2] = tmp;
  heap[i1].heap_idx = i1;
  heap[i2].heap_idx = i2;
  opt->hypernodes[heap[i1].id].heap_idx = i1;
  opt->hypernodes[heap[i2].id].heap_idx = i2;
}


void heapify(node* heap, int i, int n, optim* opt){
  if (n <= 1) return;  // No need to heapify if there's only one element
  int pos = i;
  while (2*pos + 1 < n){
    if (2*pos + 1== n-1){
      int child1_idx = 2*pos+1;
      if (heap[pos].r > heap[child1_idx].r){
        heap_swap(heap, pos, child1_idx, opt);
      }
      return;
    }
    int child1_idx = 2*pos+1;
    int child2_idx = 2*pos+2;
    if (heap[pos].r <= heap[child1_idx].r && heap[pos].r <= heap[child2_idx].r) return; 
    int min_idx = child1_idx;
    if (heap[child2_idx].r < heap[min_idx].r) min_idx = child2_idx;
    heap_swap(heap, pos, min_idx, opt);
    pos = min_idx;
  }
}



void build_heap(node* heap, int n, optim* opt){
  for (int i = n/2 -1; i >= 0; i--) heapify(heap, i, n, opt);
}

void delete_min(node* heap, int n, optim* opt){
  heap_swap(heap, 0, n-1, opt);
  heapify(heap, 0, n-1, opt);
}

void push_up(node* heap, int i, optim* opt){
  int pos = i;
  while (pos>0){
    int parent = (pos-1)/2;
    if (heap[parent].r <= heap[pos].r) return;
    heap_swap(heap, parent, pos, opt);
    pos = parent;
  }

}


//in greedy++, r is the density. In others, r is the payload. 
void onepass_Greedypp(optim *opt){
  opt->iter++;
  for (int i = 0; i < opt->n;i++){
    opt->hypernodes[i].r = opt->hypernodes[i].r + opt->hypernodes[i].d_G / opt->hypernodes[i].weight;
  }
  node* min_heap = (node*)malloc(opt->n*sizeof(node));
  int heap_size = opt->n;
  for (int i = 0; i < opt->n; i++){
     min_heap[i] = opt->hypernodes[i];
     opt->hypernodes[i].heap_idx = i;
     //min_heap[i].heap_idx = i;
     min_heap[i].r = opt->hypernodes[i].r;
     min_heap[i].heap_idx = i;
  }
  int* removed = (int*)malloc(opt->n * sizeof(int));
  for (int i=0; i < opt->n; i++) removed[i] = 0;
  int* removed_hyper_e = (int*)malloc(opt->m * sizeof(int));
  for (int i=0; i < opt->m; i++) removed_hyper_e[i] = 0;
  build_heap(min_heap, heap_size, opt);
  // for (int i = 0; i < opt->n; i++){
  //   printf("mark rate %d\n", min_heap[0].id);
  //   delete_min(min_heap, heap_size--, opt);
  // }
  // exit(-1);
//  int count = 0;
  for (int i = 0; i < opt->n; i++){
    removed[min_heap[0].id] = 1;
    //printf("min density %lf, id = %d\n", min_heap[0].r, min_heap[0].id);
    int hyper_node_id0 = min_heap[0].id;
    delete_min(min_heap, heap_size--, opt);
    // count++;
    // if (count == 1000) exit(-1);
    for (int j =0; j<opt->hypernodes[hyper_node_id0].degree;j++){
      int edge_idx = opt->hypernodes[hyper_node_id0].neighboredges[j];
      int hyperedge_idx = opt->edges[edge_idx].s;
      if (removed_hyper_e[hyperedge_idx] == 0) removed_hyper_e[hyperedge_idx] = 1;
      else continue;
      double hyperedge_weight = opt->hyperedges[hyperedge_idx].weight;
      for (int j2 = 0; j2 < opt->hyperedges[hyperedge_idx].degree; j2++){
        int edge_idx2 = opt->hyperedges[hyperedge_idx].neighboredges[j2];
        int hypernode_idx = opt->edges[edge_idx2].t;
        if (removed[hypernode_idx] ==0){
          //printf("original density0 = %lf, d_g = %lf\n", opt->hypernodes[hypernode_idx].r, opt->hypernodes[hypernode_idx].d_G );
          //printf("original density = %lf\n", min_heap[opt->hypernodes[hypernode_idx].heap_idx].r);
          //printf("ids = %d, %d\n", hypernode_idx, hyperedge_idx);
          int x1 = opt->hypernodes[hypernode_idx].heap_idx;
          //printf("density %lf\n", min_heap[x1].r);
          opt->hypernodes[hypernode_idx].r = opt->hypernodes[hypernode_idx].r - hyperedge_weight / opt->hypernodes[hypernode_idx].weight;
          min_heap[opt->hypernodes[hypernode_idx].heap_idx].r = opt->hypernodes[hypernode_idx].r;
          //printf("diff = %lf\n", hyperedge_weight / opt->hypernodes[hypernode_idx].weight);
          push_up(min_heap, opt->hypernodes[hypernode_idx].heap_idx, opt);
          //printf("current density = %lf\n", min_heap[opt->hypernodes[hypernode_idx].heap_idx].r);
          // if (opt->hypernodes[hypernode_idx].r < 0){
          //   exit(-1);
          // }
        }
      }
    }
    //printf("min id = %d\n", min_heap[0].id);
  }
  free(min_heap);
  free(removed);
  free(removed_hyper_e);
}


void onepass_Elistpp(optim *opt){

  opt->iter++;
  for (int i=0; i < opt->m; i++){
    int min_idx=0;
    int edge_idx3 = opt->hyperedges[i].neighboredges[min_idx];
    int hypernode_idx3 = opt->edges[edge_idx3].t;
    double min_r =opt->hypernodes[hypernode_idx3].r / opt->hypernodes[hypernode_idx3].weight; //min_r of the first neighbor of hypernode i. 
    for (int j=0;j<opt->hyperedges[i].degree;j++){
      int edge_idx = opt->hyperedges[i].neighboredges[j];
      int hypernode_idx = opt->edges[edge_idx].t;
      if (opt->hypernodes[hypernode_idx].r / opt->hypernodes[hypernode_idx].weight < min_r){
        min_r = opt->hypernodes[hypernode_idx].r / opt->hypernodes[hypernode_idx].weight;
        min_idx = j;
      }
    }
    int edge_idx2 = opt->hyperedges[i].neighboredges[min_idx];
    int hypernode_idx2 = opt->edges[edge_idx2].t;
    
    opt->hypernodes[hypernode_idx2].r = opt->hypernodes[hypernode_idx2].r + opt->hyperedges[i].weight;
  }
}

void Out_put_weights(optim *opt, char *file_name, int rep){
  FILE *file;
  file = fopen(file_name, "w"); // 以写入模式打开文件
  printf("%s", file_name);
  if (file == NULL) {
      printf("Cannot Open\n");
      return;
  }
  long unsigned i = 0;
  fprintf(file, "%d\n", rep);
  for (i=0; i < opt->e; i++){
    fprintf(file, "%lu, %f\n", i, opt->edges[i].a);
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

void Out_put_density_Elist(optim *opt, char *file_name, int num){
  FILE *file;
  file = fopen(file_name, "w"); // 以写入模式打开文件
  printf("%s", file_name);
  if (file == NULL) {
      printf("Cannot Open\n");
      return;
  }
  fprintf(file, "%d\n", num);
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




double multiplicative(optim *opt){
  double res= 0.0;
  for (int i=0; i < opt->n; i++){
    double error = (opt->hypernodes[i].r - opt->hypernodes[i].r_correct) / (opt->hypernodes[i].r_correct);
    if (error < 0) error = -error;
    if (error > res) res = error;
  }
  return res;
}

double absolute(optim *opt){
  double res = 0.0;
  for (int i=0; i < opt->n; i++){
    double error = (opt->hypernodes[i].r - opt->hypernodes[i].r_correct) * (opt->hypernodes[i].r - opt->hypernodes[i].r_correct) / (opt->hypernodes[i].weight * opt->hypernodes[i].weight);
    res = res + error;
  }
  res = sqrt(res);
  return res;

}

double absolute2(optim *opt, int T){
  double res = 0.0;
  for (int i=0; i < opt->n; i++){
    double error = (opt->hypernodes[i].r / T - opt->hypernodes[i].r_correct) * (opt->hypernodes[i].r /T - opt->hypernodes[i].r_correct) / (opt->hypernodes[i].weight * opt->hypernodes[i].weight);
    res = res + error;
  }
  res = sqrt(res);
  return res;
}

double absolute_greedy(optim *opt, int T){
  double res = 0.0;
  for (int i=0; i < opt->n; i++){
    double error = (opt->hypernodes[i].r * opt->hypernodes[i].weight / T - opt->hypernodes[i].r_correct) * (opt->hypernodes[i].r * opt->hypernodes[i].weight /T - opt->hypernodes[i].r_correct) / (opt->hypernodes[i].weight * opt->hypernodes[i].weight);
    res = res + error;
  }
  res = sqrt(res);
  return res;
}

double multiplicative2(optim *opt, int T){
  double res= 0.0;
  for (int i=0; i < opt->n; i++){
    double error = (opt->hypernodes[i].r / T - opt->hypernodes[i].r_correct) / (opt->hypernodes[i].r_correct);
    if (error < 0) error = -error;
    if (error > res) res = error;
  }
  return res;
}

double multiplicative_greedy(optim *opt, int T){
  double res= 0.0;
  for (int i=0; i < opt->n; i++){
    double error = (opt->hypernodes[i].r * opt->hypernodes[i].weight / T - opt->hypernodes[i].r_correct) / (opt->hypernodes[i].r_correct);
    if (error < 0) error = -error;
    if (error > res) res = error;
  }
  return res;
}

long long merge(node* Mark_nodes, int l, int mid, int r){
  long long res=0;
  node* temp_arr = (node*)malloc((mid-l+1)*sizeof(node));
  for (int i = 0; i < mid - l + 1; i++) temp_arr[i].r_correct = Mark_nodes[i+l].r_correct;
  int i1=0;
  int i2 = mid+1;
  int j=l;
  while (i1 <mid - l + 1 && i2<=r){
    if (temp_arr[i1].r_correct <= Mark_nodes[i2].r_correct){
      Mark_nodes[j++] = temp_arr[i1++];
    }
    else{
      res = res + (mid-l - i1 + 1);
      Mark_nodes[j++] = Mark_nodes[i2++];
    }
  }
  while (i1 < mid - l + 1) Mark_nodes[j++] = temp_arr[i1++];
  free(temp_arr);
  return res;
}


long long num_inversions(node* Mark_nodes, int l, int r){
  if (l >= r) return 0;
  long long mid = (l+r)/2;
  long long res1 = num_inversions(Mark_nodes, l, mid);
  long long res2 = num_inversions(Mark_nodes, mid+1, r);
  long long res3 = merge(Mark_nodes, l, mid, r);
  return res1+res2+res3;

}

//used for quicksort, smaller r has higher priority. Be careful here.  
static int compare_nodes2(void const *a, void const *b){
	node const *pa = (const node*)a;
	node const *pb = (const node*)b;
	if ((*pa).r<(*pb).r)
		return -1;
  else if ((*pa).r>(*pb).r){
    return 1;
  }
  else if ((*pa).id<(*pb).id){
    return -1;
  }
  else return 1;
	return -1;
}


void initialize(optim *opt){
  opt->iter = 0;
	for (int k=0;k<opt->e;k++){
    int hyperedge_idx = opt->edges[k].s;
    double weight = opt->hyperedges[hyperedge_idx].weight;
    int deg = opt->hyperedges[hyperedge_idx].degree;
    opt->edges[k].a= weight / deg;
  }
  update_density(opt);
}



// 获取当前时间（毫秒）
long long getCurrentTimeInMilliseconds() {
    struct timeval currentTime;
    gettimeofday(&currentTime, NULL);
    return currentTime.tv_sec * 1000LL + currentTime.tv_usec / 1000;
}

double obj_function(optim *opt){
  update_density(opt);
  double res = 0.0;
  for (int i = 0; i < opt->n; i++){
    res = res + opt->hypernodes[i].r * opt->hypernodes[i].r / opt->hypernodes[i].weight;
    
  }
  return res;
}

double obj_function_greedy(optim *opt, int T){
  //update_density(opt);
  double res = 0.0;
  for (int i = 0; i < opt->n; i++){
    res = res + opt->hypernodes[i].r / T * opt->hypernodes[i].r / T * opt->hypernodes[i].weight;
  }
  return res;
}


double obj_function2(optim *opt, int T){
  //update_density(opt);
  double res = 0.0;
  for (int i = 0; i < opt->n; i++){
    res = res + opt->hypernodes[i].r / T * opt->hypernodes[i].r / T / opt->hypernodes[i].weight;
  }
  return res;
}

void read_weights(optim *opt, char *name1, char *name2){
  //------read edge part-----------------------
  FILE *file;
  
	file=fopen(name1,"r");
  int idx;
  double weight;
	while (fscanf(file,"%d %lf", &(idx), &(weight))==2) {
    //opt->edges[opt->e].s = 
		int new_idx = opt->newlabel[idx];
    opt->hyperedges[new_idx].weight = weight;
  }
	fclose(file);

	file=fopen(name2,"r");
	while (fscanf(file,"%d %lf", &(idx), &(weight))==2) {
    //opt->edges[opt->e].s = 
		int new_idx = opt->newlabel2[idx];
    opt->hypernodes[new_idx].weight = weight;
  }
	fclose(file);
  for (int k=0;k<opt->n;k++){
    opt->hypernodes[k].d_G = 0.0;
  }
  for (int k=0;k<opt->n;k++){
    for (int x=0; x < opt->hypernodes[k].degree; x++){
      int edge_idx = opt->hypernodes[k].neighboredges[x];
      int hyperedge_idx = opt->edges[edge_idx].s;
      double weight = opt->hyperedges[hyperedge_idx].weight;
      opt->hypernodes[k].d_G = opt->hypernodes[k].d_G + weight;
    }
    // if (opt->hypernodes[k].d_G - opt->hypernodes[k].degree > 0.00000001 || opt->hypernodes[k].d_G - opt->hypernodes[k].degree < -0.00000001){
    //   printf("Strange!!!\n");
    // }
  }
}
using namespace std;
int main(int argc,char** argv){
	optim* opt;
  //isoreg *fit,*fit2;
	unsigned nsgs,i,k;
	unsigned nthreads=4;
	unsigned rep=1000;
   
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
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Building the datastructure\n");
	t1=time(NULL);
	relabel(opt); 

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	printf("- Building the datastructure\n");
	printf("- Number of nodes = %u\n",opt->n);
	printf("- Number of edges = %llu\n",opt->e);
	printf("- Computing the locally densest decomposition\n");
  printf("- Step 1: Frank-Wolf gradiant descent (%u iterations)\n",rep);
	t1=time(NULL);
	init(opt);
  //read_weights(opt, "bipartite_mark_E.txt", "bipartite_mark_N.txt");
  // if we use FISTA, we need the following. Otherwise, we can comment them. 
  double max_degree = get_max_weighted_degree(opt);
	link0 *z; 
	z = (link0*)malloc(opt->e*sizeof(link0));
	link0 *y; 
	y = (link0*)malloc(opt->e*sizeof(link0));
	link0 *last_x;
	last_x = (link0*)malloc(opt->e*sizeof(link0));
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
  //-----file name====================
  double real_norm = 0.0; //-------store correct norm...
  double real_obj = 0.0; //------------store correct obj function...
  printf("Read the correct answer\n");
  FILE *exact = fopen("file_density.txt", "r");
  int id = 0;
  double rho = 0.0;
	while (fscanf(exact,"%d, %lf", &(id), &(rho))==2) {
    opt->hypernodes[id].r_correct = rho * opt->hypernodes[id].weight;
    real_norm = real_norm + opt->hypernodes[id].r_correct * opt->hypernodes[id].r_correct;
    real_obj = real_obj + opt->hypernodes[id].r_correct * opt->hypernodes[id].r_correct / opt->hypernodes[id].weight;
    //printf("density is %lf\n", rho);
	}
  real_norm = sqrt(real_norm);
    fclose(exact);
  printf("real norm%lf\n", real_norm);
  node *Mark_nodes = (node*)malloc(opt->n * sizeof(node));
  printf("start the algorithm...\n");
  //Out_put_density(opt, "output_init.txt");
  //---------------------------------------
  double lambda = 0;
  int flag_pp = 1;
  int num_points = 1;
  long int total_time = 0;
  double average_time = 0.0;
  double** simulation = (double**)malloc(rep*sizeof(double*));
  long long* inv_arr = (long long*)malloc(rep*sizeof(long long));
  for (i=0; i < rep; i++) simulation[i] = (double*)malloc(7*sizeof(double));
  for (i=0; i < rep; i++) simulation[i] = (double*)malloc(7*sizeof(double));
  FILE* file_error;
  FILE* file_time;
  FILE* file_abs_error;
  FILE* file_abs_time;
  FILE* file_abs_normalized;
  FILE* file_abs_normalized_time;

  FILE* file_inv;
  FILE* file_inv_time;
  
  
  FILE* file_obj;
  FILE* file_obj_normalized;
  FILE* file_obj_time;
  FILE* file_obj_normalized_time;

  FILE* file_Abserror;
  //int num = 0;
  double error = 0;
  double abs_error = 0;
  long long inv = 0;
  double obj = 0;
  int num = 0;

  long long int time1 = 0;
  long long int time2 = 0;
  string mul_name = "_mul.txt";
  string mul_time_name = "_mul_time.txt";
  string inv_name = "_inv.txt";
  string inv_time_name = "_inv_time.txt";
  string abs_name = "_abs.txt";
  string abs_time_name = "_abs_time.txt";
  string abs_normalized_name = "_abs_normalized.txt";
  string abs_normalized_time_name = "_abs_normalized_time.txt";
  string obj_name = "_obj.txt";
  string obj_time_name = "_obj_time.txt";
  string obj_norm_name = "_obj_normalized.txt";
  string obj_norm_time_name = "_obj_normalized_time.txt";
  string algo_name = "";
  string folder_name = "data_normal/";
  string file_name_cpp = "";
  const char* file_name_c;
  int flag_increase = 0;
  vector<string> algos = {}; //"PR_exp", "FISTA*"
  // Iterate over the vector of algorithms "FISTA", , "PR_exp" "FISTA*", 
  for (const std::string& algo_name : algos) {
    initialize(opt);

    // 如果是PR+ momentum, 计算momentum 所需要的东西。
    int flag_init = 0;
    size_t found = algo_name.find("lin");
    if (found != string::npos) flag_init = 1;
    found = algo_name.find("exp");
    if (found != string::npos) flag_init = 1;
    if (flag_init == 1){
      //printf("mark calculate\n");
      for (int i1 = 0; i1 < opt->e; i1++){
        last_x[i1].a = opt->edges[i1].a;
        z[i1].a = opt->edges[i1].a;
        //last_x[i].b = 1 - opt->edges[i].a;
      }
    }
    //如果是fista，计算。
    if (algo_name == "FISTA" || algo_name == "FISTA*"){
      for (int i1 = 0; i1 < opt->e; i1++){
        last_x[i1].a = opt->edges[i1].a;
        z[i1].a = opt->edges[i1].a;
        y[i1].a = opt->edges[i1].a;
        //last_x[i].b = 1 - opt->edges[i].a;
      }
    }
    
    //algo_name = "PR_exp";
    total_time = 0;
    cout << "start " +algo_name  + " algorithm..." << endl;
    //multiplicative error and time--------------------------
    file_name_cpp = folder_name + algo_name + mul_name;
    file_name_c = file_name_cpp.c_str();
    file_error = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + mul_time_name;
    file_name_c = file_name_cpp.c_str();
    file_time = fopen(file_name_c, "a");
    //absolute error and time--------------------------
    file_name_cpp = folder_name + algo_name + abs_name;
    file_name_c = file_name_cpp.c_str();
    file_abs_error = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + abs_time_name;
    file_name_c = file_name_cpp.c_str();
    file_abs_time = fopen(file_name_c, "a");
    //normalized absolute error and time--------------------------
    file_name_cpp = folder_name + algo_name + abs_normalized_name;
    file_name_c = file_name_cpp.c_str();
    file_abs_normalized = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
    file_name_c = file_name_cpp.c_str();
    file_abs_normalized_time = fopen(file_name_c, "a");
    //number of inversions and time-------------------------
    file_name_cpp = folder_name + algo_name + inv_name;
    file_name_c = file_name_cpp.c_str();
    file_inv = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + inv_time_name;
    file_name_c = file_name_cpp.c_str();
    file_inv_time = fopen(file_name_c, "a");
    //objective function and time/ normalized objtive function and time.-------------------------
    file_name_cpp = folder_name + algo_name + obj_name;
    file_name_c = file_name_cpp.c_str();
    file_obj = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + obj_norm_name;
    file_name_c = file_name_cpp.c_str();
    file_obj_normalized = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + obj_time_name;
    file_name_c = file_name_cpp.c_str();
    file_obj_time = fopen(file_name_c, "a");
    file_name_cpp = folder_name + algo_name + obj_norm_time_name;
    file_name_c = file_name_cpp.c_str();
    file_obj_normalized_time = fopen(file_name_c, "a");
    
    // double error = 0;
    // struct timeval Time1;
    // struct timeval Time2;


    error = multiplicative(opt);
    abs_error = absolute(opt);
    obj = obj_function(opt);

    //printf("the error is %lf\n", error);
    simulation[0][0] = error; //multiplicative error
    simulation[0][1] = (double)(0);
    simulation[0][2] = abs_error; //absolute error
    simulation[0][3] = abs_error / real_norm;
    simulation[0][4] = obj;
    simulation[0][5] = obj - real_obj;
    fprintf(file_error, "%d %lf\n", (0), error);
    fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
    fprintf(file_abs_error, "%d %lf\n", (0), simulation[0][2]);
    for (int i1 =0; i1 < opt->n; i1++){
      Mark_nodes[i1] = opt->hypernodes[i1];
    }
    qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
    inv = num_inversions(Mark_nodes, 0, opt->n-1);
    inv_arr[0] = inv;
    fprintf(file_inv, "%d %lld\n", (0), inv);
    fprintf(file_obj_normalized, "%d %lf\n", (0), simulation[0][5]);
    fprintf(file_abs_normalized, "%d %lf\n", (0), simulation[0][3]);
    //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
    
    for (i=0;i<rep;i++){
      //printf("%u\n",i);
      //onepass(opt);
      //onepass_PR_plus(opt, last_x, z, &lambda);
      //onepass_PR(opt);
      //onepass_Greedypp(opt);
      //onepass_Elistpp(opt);
      time1 = getCurrentTimeInMilliseconds();
      printf("start iteration %d\n", i);
      //onepass_PR_plus_normal(opt, last_x, z, &lambda);
      //onepass_PR_plus_new_2_normal_change(opt, last_x, z, last_gradient);
      //{"FW_WF", "PR_exp", "FW_syn", "FISTA", "FW_Elist", "PR", "PR_lin", "PR_exp_pro"};
      if (algo_name == "FW_WF") onepass_FW_Water(opt);
      else if (algo_name == "PR_exp_lin") onepass_PR_plus_product_normal(opt, last_x, z);
      else if (algo_name == "FW_syn") onepass_FW_syn(opt);
      else if (algo_name == "FISTA") onepass_FISTA_normal(opt, max_degree, z, y, last_x);
      else if (algo_name == "FW_Elist") onepass_FW_Elist(opt);
      else if (algo_name == "PR") onepass_PR(opt);
      else if (algo_name == "PR_lin") onepass_PR_plus_lin_normal(opt, last_x, z);
      else if (algo_name == "PR_exp") onepass_PR_plus_product_normal_pro(opt, last_x, z);
      else if (algo_name == "FISTA*") onepass_FISTA_star_normal(opt, max_degree, z, y, last_x);
      //else if ((algo_name == "FISTA*") onepass_FISTA_star_normal(opt, last_x, z);)
      //onepass_PR_plus_product_normal(opt, last_x, z);
      //onepass_PR_plus_new_2(opt, last_x, z);
      time2 = getCurrentTimeInMilliseconds();
        //test
      // for (int i1 = 0; i1 < opt->n; i1++){
      //   printf("r = %lf\n", opt->hypernodes[i1].r);
      // }
      // scanf("%d", &rep);
      total_time = total_time + time2 - time1;
      //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
      //onepass_FISTA(opt, max_degree, z, y, last_x);
      //printf("iteration %d finished\n", i+1);
      //printf("the error is %lf\n", error);
      if ((i+1)%num_points == 0){
        error = multiplicative(opt);
        abs_error = absolute(opt);
        obj = obj_function(opt);

        //printf("the error is %lf\n", error);
        simulation[i][0] = error; //multiplicative error
        simulation[i][1] = (double)(i+1);
        simulation[i][2] = abs_error; //absolute error
        simulation[i][3] = abs_error / real_norm;
        simulation[i][4] = obj;
        simulation[i][5] = obj - real_obj;
        fprintf(file_error, "%d %lf\n", (i+1), error);
        fprintf(file_obj, "%d %lf\n", (i+1), obj);
        fprintf(file_abs_error, "%d %lf\n", (i+1), simulation[i][2]);
        for (int i1 =0; i1 < opt->n; i1++){
          Mark_nodes[i1] = opt->hypernodes[i1];
        }
        qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
        inv = num_inversions(Mark_nodes, 0, opt->n-1);
        inv_arr[i] = inv;
        if (inv == 0) flag_increase = 0;
        fprintf(file_inv, "%d %lld\n", (i+1), inv);
        fprintf(file_obj_normalized, "%d %lf\n", (i+1), simulation[i][5]);
        fprintf(file_abs_normalized, "%d %lf\n", (i+1), simulation[i][3]);
      }
      //output intermediate weights. 
      if (i % 100 == 0){
        file_name_cpp = algo_name+"weight.txt";
        file_name_c = file_name_cpp.c_str();
        Out_put_weights(opt, (char *)file_name_c, i);
        fclose(file_error);
        fclose(file_time);
        fclose(file_abs_error);
        fclose(file_abs_time);
        fclose(file_abs_normalized);
        fclose(file_abs_normalized_time);
        fclose(file_obj);
        fclose(file_obj_time);
        fclose(file_obj_normalized);
        fclose(file_obj_normalized_time);
        fclose(file_inv);
        fclose(file_inv_time);
        //multiplicative error and time--------------------------
        file_name_cpp = folder_name + algo_name + mul_name;
        file_name_c = file_name_cpp.c_str();
        file_error = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + mul_time_name;
        file_name_c = file_name_cpp.c_str();
        file_time = fopen(file_name_c, "a");
        //absolute error and time--------------------------
        file_name_cpp = folder_name + algo_name + abs_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_error = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + abs_time_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_time = fopen(file_name_c, "a");
        //normalized absolute error and time--------------------------
        file_name_cpp = folder_name + algo_name + abs_normalized_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_normalized = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_normalized_time = fopen(file_name_c, "a");
        //number of inversions and time-------------------------
        file_name_cpp = folder_name + algo_name + inv_name;
        file_name_c = file_name_cpp.c_str();
        file_inv = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + inv_time_name;
        file_name_c = file_name_cpp.c_str();
        file_inv_time = fopen(file_name_c, "a");
        //objective function and time/ normalized objtive function and time.-------------------------
        file_name_cpp = folder_name + algo_name + obj_name;
        file_name_c = file_name_cpp.c_str();
        file_obj = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + obj_norm_name;
        file_name_c = file_name_cpp.c_str();
        file_obj_normalized = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + obj_time_name;
        file_name_c = file_name_cpp.c_str();
        file_obj_time = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + obj_norm_time_name;
        file_name_c = file_name_cpp.c_str();
        file_obj_normalized_time = fopen(file_name_c, "a");
        }
        // if (i>=rep-1 && flag_increase == 1){
        //   rep = rep + 100;
        //   simulation = (double**)realloc(simulation, rep*sizeof(double*));
        //   inv_arr = (long long*)realloc(inv_arr, rep*sizeof(long long));
        //   for (int i0=rep-100; i0 < rep; i0++) simulation[i0] = (double*)malloc(7*sizeof(double));
        //   Out_put_weights(opt, "weights.txt", i);
        // }
	  
    }


    average_time = (double)total_time / (double)rep;
    average_time = average_time / 1000;
    fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
    fprintf(file_abs_time, "%lf %lf\n", simulation[0][2], 0.0);
    fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[0][3], 0.0);
    fprintf(file_obj_time, "%lf %lf\n", simulation[0][4], 0.0);
    fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[0][5], 0.0);
    fprintf(file_inv_time, "%lld %lf\n", inv_arr[0], 0.0);

    for (i=0;i<rep;i++){
      if ((i+1)%num_points == 0){
        double sim_time = simulation[i][1]*average_time;
        fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
        fprintf(file_abs_time, "%lf %lf\n", simulation[i][2], sim_time);
        fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[i][3], sim_time);
        fprintf(file_obj_time, "%lf %lf\n", simulation[i][4], sim_time);
        fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[i][5], sim_time);
        fprintf(file_inv_time, "%lld %lf\n", inv_arr[i], sim_time);
        if (sim_time < 0){
          printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
        }
      }
    }



    fclose(file_error);
    fclose(file_time);
    fclose(file_abs_error);
    fclose(file_abs_time);
    fclose(file_abs_normalized);
    fclose(file_abs_normalized_time);
    fclose(file_obj);
    fclose(file_obj_time);
    fclose(file_obj_normalized);
    fclose(file_obj_normalized_time);
    fclose(file_inv);
    fclose(file_inv_time);
    
  }



  initialize(opt);
  for (int i = 0; i < opt->n;i++){
    opt->hypernodes[i].r = 0.0;
  }
  total_time = 0;
  algo_name = "Elist++";
  cout << "start " +algo_name  + " algorithm..." << endl;
  //multiplicative error and time--------------------------
  file_name_cpp = folder_name + algo_name + mul_name;
  file_name_c = file_name_cpp.c_str();
  file_error = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + mul_time_name;
  file_name_c = file_name_cpp.c_str();
  file_time = fopen(file_name_c, "w");
  //absolute error and time--------------------------
  file_name_cpp = folder_name + algo_name + abs_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_error = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + abs_time_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_time = fopen(file_name_c, "w");
  //normalized absolute error and time--------------------------
  file_name_cpp = folder_name + algo_name + abs_normalized_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_normalized = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_normalized_time = fopen(file_name_c, "w");
  //number of inversions and time-------------------------
  file_name_cpp = folder_name + algo_name + inv_name;
  file_name_c = file_name_cpp.c_str();
  file_inv = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + inv_time_name;
  file_name_c = file_name_cpp.c_str();
  file_inv_time = fopen(file_name_c, "w");
  //objective function and time/ normalized objtive function and time.-------------------------
  file_name_cpp = folder_name + algo_name + obj_name;
  file_name_c = file_name_cpp.c_str();
  file_obj = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + obj_norm_name;
  file_name_c = file_name_cpp.c_str();
  file_obj_normalized = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + obj_time_name;
  file_name_c = file_name_cpp.c_str();
  file_obj_time = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + obj_norm_time_name;
  file_name_c = file_name_cpp.c_str();
  file_obj_normalized_time = fopen(file_name_c, "w");
  
  // double error = 0;
  // struct timeval Time1;
  // struct timeval Time2;


  //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, z, &lambda);
		//onepass_PR(opt);
    //onepass_Greedypp(opt);
    //onepass_Elistpp(opt);
    time1 = getCurrentTimeInMilliseconds();
    //onepass_PR_plus_normal(opt, last_x, z, &lambda);
		//onepass_PR_plus_new_2_normal_change(opt, last_x, z, last_gradient);
    //onepass_Greedypp(opt);
    onepass_Elistpp(opt);
    //onepass_PR_plus_new_2(opt, last_x, z);
    time2 = getCurrentTimeInMilliseconds();
      //test
    // for (int i1 = 0; i1 < opt->n; i1++){
    //   printf("r = %lf\n", opt->hypernodes[i1].r);
    // }
    // scanf("%d", &rep);
    total_time = total_time + time2 - time1;
    //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
    //printf("iteration %d finished\n", i+1);
    //printf("the error is %lf\n", error);
    if ((i+1)%num_points == 0){
      error = multiplicative2(opt, i+1);
      abs_error = absolute2(opt, i+1);
      obj = obj_function2(opt, i+1);

      //printf("the error is %lf\n", error);
      simulation[i][0] = error; //multiplicative error
      simulation[i][1] = (double)(i+1);
      simulation[i][2] = abs_error; //absolute error
      simulation[i][3] = abs_error / real_norm;
      simulation[i][4] = obj;
      simulation[i][5] = obj - real_obj;
      fprintf(file_error, "%d %lf\n", (i+1), error);
      fprintf(file_obj, "%d %lf\n", (i+1), obj);
      fprintf(file_abs_error, "%d %lf\n", (i+1), simulation[i][2]);
      for (int i1 =0; i1 < opt->n; i1++){
        Mark_nodes[i1] = opt->hypernodes[i1];
      }
      qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
      inv = num_inversions(Mark_nodes, 0, opt->n-1);
      inv_arr[i] = inv;
      fprintf(file_inv, "%d %lld\n", (i+1), inv);
      fprintf(file_obj_normalized, "%d %lf\n", (i+1), simulation[i][5]);
      fprintf(file_abs_normalized, "%d %lf\n", (i+1), simulation[i][3]);
    }

    if ((i+1) % 1 == 0){
        //Out_put_density_Elist(opt, "Elist++.txt", i);
        fclose(file_error);
        fclose(file_time);
        fclose(file_abs_error);
        fclose(file_abs_time);
        fclose(file_abs_normalized);
        fclose(file_abs_normalized_time);
        fclose(file_obj);
        fclose(file_obj_time);
        fclose(file_obj_normalized);
        fclose(file_obj_normalized_time);
        fclose(file_inv);
        fclose(file_inv_time);
        exit(-1);
        //multiplicative error and time--------------------------
        file_name_cpp = folder_name + algo_name + mul_name;
        file_name_c = file_name_cpp.c_str();
        file_error = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + mul_time_name;
        file_name_c = file_name_cpp.c_str();
        file_time = fopen(file_name_c, "a");
        //absolute error and time--------------------------
        file_name_cpp = folder_name + algo_name + abs_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_error = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + abs_time_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_time = fopen(file_name_c, "a");
        //normalized absolute error and time--------------------------
        file_name_cpp = folder_name + algo_name + abs_normalized_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_normalized = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
        file_name_c = file_name_cpp.c_str();
        file_abs_normalized_time = fopen(file_name_c, "a");
        //number of inversions and time-------------------------
        file_name_cpp = folder_name + algo_name + inv_name;
        file_name_c = file_name_cpp.c_str();
        file_inv = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + inv_time_name;
        file_name_c = file_name_cpp.c_str();
        file_inv_time = fopen(file_name_c, "a");
        //objective function and time/ normalized objtive function and time.-------------------------
        file_name_cpp = folder_name + algo_name + obj_name;
        file_name_c = file_name_cpp.c_str();
        file_obj = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + obj_norm_name;
        file_name_c = file_name_cpp.c_str();
        file_obj_normalized = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + obj_time_name;
        file_name_c = file_name_cpp.c_str();
        file_obj_time = fopen(file_name_c, "a");
        file_name_cpp = folder_name + algo_name + obj_norm_time_name;
        file_name_c = file_name_cpp.c_str();
        file_obj_normalized_time = fopen(file_name_c, "a");
    }


	}


  average_time = (double)total_time / (double)rep;
  average_time = average_time / 1000;
  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_abs_time, "%lf %lf\n", simulation[0][2], 0.0);
  // fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[0][3], 0.0);
  // fprintf(file_obj_time, "%lf %lf\n", simulation[0][4], 0.0);
  // fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[0][5], 0.0);
  // fprintf(file_inv_time, "%lf %lf\n", simulation[0][6], 0.0);

	for (i=0;i<rep;i++){
    if ((i+1)%num_points == 0){
      double sim_time = simulation[i][1]*average_time;
      fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
      fprintf(file_abs_time, "%lf %lf\n", simulation[i][2], sim_time);
      fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[i][3], sim_time);
      fprintf(file_obj_time, "%lf %lf\n", simulation[i][4], sim_time);
      fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[i][5], sim_time);
      fprintf(file_inv_time, "%lld %lf\n", inv_arr[i], sim_time);
      if (sim_time < 0){
        printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
      }
    }
	}



  fclose(file_error);
  fclose(file_time);
  fclose(file_abs_error);
  fclose(file_abs_time);
  fclose(file_abs_normalized);
  fclose(file_abs_normalized_time);
  fclose(file_obj);
  fclose(file_obj_time);
  fclose(file_obj_normalized);
  fclose(file_obj_normalized_time);
  fclose(file_inv);
  fclose(file_inv_time);



  return 1;


  // for (int i=0; i < opt->n; i++){
  //   printf("rate = %lf\n", Mark_nodes[i].r_correct);
  // }

}


/*


//--------------------specical algorithms---------------------------------------------------

  initialize(opt);
  for (int i = 0; i < opt->n;i++){
    opt->hypernodes[i].r = 0.0;
  }
  total_time = 0;
  algo_name = "Greedy++";
  cout << "start " +algo_name  + " algorithm..." << endl;
  //multiplicative error and time--------------------------
  file_name_cpp = folder_name + algo_name + mul_name;
  file_name_c = file_name_cpp.c_str();
  file_error = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + mul_time_name;
  file_name_c = file_name_cpp.c_str();
  file_time = fopen(file_name_c, "w");
  //absolute error and time--------------------------
  file_name_cpp = folder_name + algo_name + abs_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_error = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + abs_time_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_time = fopen(file_name_c, "w");
  //normalized absolute error and time--------------------------
  file_name_cpp = folder_name + algo_name + abs_normalized_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_normalized = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
  file_name_c = file_name_cpp.c_str();
  file_abs_normalized_time = fopen(file_name_c, "w");
  //number of inversions and time-------------------------
  file_name_cpp = folder_name + algo_name + inv_name;
  file_name_c = file_name_cpp.c_str();
  file_inv = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + inv_time_name;
  file_name_c = file_name_cpp.c_str();
  file_inv_time = fopen(file_name_c, "w");
  //objective function and time/ normalized objtive function and time.-------------------------
  file_name_cpp = folder_name + algo_name + obj_name;
  file_name_c = file_name_cpp.c_str();
  file_obj = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + obj_norm_name;
  file_name_c = file_name_cpp.c_str();
  file_obj_normalized = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + obj_time_name;
  file_name_c = file_name_cpp.c_str();
  file_obj_time = fopen(file_name_c, "w");
  file_name_cpp = folder_name + algo_name + obj_norm_time_name;
  file_name_c = file_name_cpp.c_str();
  file_obj_normalized_time = fopen(file_name_c, "w");
  
  // double error = 0;
  // struct timeval Time1;
  // struct timeval Time2;


  //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, z, &lambda);
		//onepass_PR(opt);
    //onepass_Greedypp(opt);
    //onepass_Elistpp(opt);
    time1 = getCurrentTimeInMilliseconds();
    //onepass_PR_plus_normal(opt, last_x, z, &lambda);
		//onepass_PR_plus_new_2_normal_change(opt, last_x, z, last_gradient);
    onepass_Greedypp(opt);
    //onepass_PR_plus_new_2(opt, last_x, z);
    time2 = getCurrentTimeInMilliseconds();
      //test
    // for (int i1 = 0; i1 < opt->n; i1++){
    //   printf("r = %lf\n", opt->hypernodes[i1].r);
    // }
    // scanf("%d", &rep);
    total_time = total_time + time2 - time1;
    //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
    //printf("iteration %d finished\n", i+1);
    //printf("the error is %lf\n", error);
    if ((i+1)%num_points == 0){
      error = multiplicative_greedy(opt, i+1);
      abs_error = absolute_greedy(opt, i+1);
      obj = obj_function_greedy(opt, i+1);

      //printf("the error is %lf\n", error);
      simulation[i][0] = error; //multiplicative error
      simulation[i][1] = (double)(i+1);
      simulation[i][2] = abs_error; //absolute error
      simulation[i][3] = abs_error / real_norm;
      simulation[i][4] = obj;
      simulation[i][5] = obj - real_obj;
      fprintf(file_error, "%d %lf\n", (i+1), error);
      fprintf(file_obj, "%d %lf\n", (i+1), obj);
      fprintf(file_abs_error, "%d %lf\n", (i+1), simulation[i][2]);
      for (int i1 =0; i1 < opt->n; i1++){
        Mark_nodes[i1] = opt->hypernodes[i1];
        Mark_nodes[i1].r = opt->hypernodes[i1].r * Mark_nodes[i1].weight;
      }
      qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
      inv = num_inversions(Mark_nodes, 0, opt->n-1);
      inv_arr[i] = inv;
      fprintf(file_inv, "%d %lld\n", (i+1), inv);
      fprintf(file_obj_normalized, "%d %lf\n", (i+1), simulation[i][5]);
      fprintf(file_abs_normalized, "%d %lf\n", (i+1), simulation[i][3]);
    }
	}


  average_time = (double)total_time / (double)rep;
  average_time = average_time / 1000;
  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_abs_time, "%lf %lf\n", simulation[0][2], 0.0);
  // fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[0][3], 0.0);
  // fprintf(file_obj_time, "%lf %lf\n", simulation[0][4], 0.0);
  // fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[0][5], 0.0);
  // fprintf(file_inv_time, "%lf %lf\n", simulation[0][6], 0.0);

	for (i=0;i<rep;i++){
    if ((i+1)%num_points == 0){
      double sim_time = simulation[i][1]*average_time;
      fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
      fprintf(file_abs_time, "%lf %lf\n", simulation[i][2], sim_time);
      fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[i][3], sim_time);
      fprintf(file_obj_time, "%lf %lf\n", simulation[i][4], sim_time);
      fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[i][5], sim_time);
      fprintf(file_inv_time, "%lld %lf\n", inv_arr[i], sim_time);
      if (sim_time < 0){
        printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
      }
    }
	}



  fclose(file_error);
  fclose(file_time);
  fclose(file_abs_error);
  fclose(file_abs_time);
  fclose(file_abs_normalized);
  fclose(file_abs_normalized_time);
  fclose(file_obj);
  fclose(file_obj_time);
  fclose(file_obj_normalized);
  fclose(file_obj_normalized_time);
  fclose(file_inv);
  fclose(file_inv_time);


*/


  // initialize(opt);
  // total_time = 0;
  // algo_name = "FW_WF";
  // cout << "start " +algo_name  + " algorithm..." << endl;
  // //multiplicative error and time--------------------------
  // file_name_cpp = folder_name + algo_name + mul_name;
  // file_name_c = file_name_cpp.c_str();
  // file_error = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + mul_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_time = fopen(file_name_c, "w");
  // //absolute error and time--------------------------
  // file_name_cpp = folder_name + algo_name + abs_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_error = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + abs_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_time = fopen(file_name_c, "w");
  // //normalized absolute error and time--------------------------
  // file_name_cpp = folder_name + algo_name + abs_normalized_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_normalized = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_normalized_time = fopen(file_name_c, "w");
  // //number of inversions and time-------------------------
  // file_name_cpp = folder_name + algo_name + inv_name;
  // file_name_c = file_name_cpp.c_str();
  // file_inv = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + inv_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_inv_time = fopen(file_name_c, "w");
  // //objective function and time/ normalized objtive function and time.-------------------------
  // file_name_cpp = folder_name + algo_name + obj_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_norm_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_normalized = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_time = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_norm_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_normalized_time = fopen(file_name_c, "w");
  
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;



  // error = multiplicative(opt);
  // abs_error = absolute(opt);
  // obj = obj_function(opt);

  // //printf("the error is %lf\n", error);
  // simulation[0][0] = error; //multiplicative error
  // simulation[0][1] = (double)(0);
  // simulation[0][2] = abs_error; //absolute error
  // simulation[0][3] = abs_error / real_norm;
  // simulation[0][4] = obj;
  // simulation[0][5] = obj - real_obj;
  // fprintf(file_error, "%d %lf\n", (0), error);
  // fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
  // fprintf(file_abs_error, "%d %lf\n", (0), simulation[0][2]);
  // for (int i1 =0; i1 < opt->n; i1++){
  //   Mark_nodes[i1] = opt->hypernodes[i1];
  // }
  // qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  // inv = num_inversions(Mark_nodes, 0, opt->n-1);
  // simulation[0][6] = inv;
  // fprintf(file_inv, "%d %d\n", (0), inv);
  // fprintf(file_obj_normalized, "%d %lf\n", (0), simulation[0][5]);
  // fprintf(file_abs_normalized, "%d %lf\n", (0), simulation[0][3]);
  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
  //   //onepass_PR_plus_normal(opt, last_x, z, &lambda);
	// 	//onepass_PR_plus_new_2_normal_change(opt, last_x, z, last_gradient);
  //   onepass_FW_Water(opt);
  //   //onepass_PR_plus_new_2(opt, last_x, z);
  //   time2 = getCurrentTimeInMilliseconds();
  //     //test
  //   // for (int i1 = 0; i1 < opt->n; i1++){
  //   //   printf("r = %lf\n", opt->hypernodes[i1].r);
  //   // }
  //   // scanf("%d", &rep);
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative(opt);
  //     abs_error = absolute(opt);
  //     obj = obj_function(opt);

  //     //printf("the error is %lf\n", error);
  //     simulation[i][0] = error; //multiplicative error
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = abs_error; //absolute error
  //     simulation[i][3] = abs_error / real_norm;
  //     simulation[i][4] = obj;
  //     simulation[i][5] = obj - real_obj;
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj);
  //     fprintf(file_abs_error, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     inv = num_inversions(Mark_nodes, 0, opt->n-1);
  //     simulation[i][6] = inv;
  //     fprintf(file_inv, "%d %d\n", (i+1), inv);
  //     fprintf(file_obj_normalized, "%d %lf\n", (i+1), simulation[i][5]);
  //     fprintf(file_abs_normalized, "%d %lf\n", (i+1), simulation[i][3]);
      
  //   }
	// }


  // average_time = (double)total_time / (double)rep;
  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_abs_time, "%lf %lf\n", simulation[0][2], 0.0);
  // fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[0][3], 0.0);
  // fprintf(file_obj_time, "%lf %lf\n", simulation[0][4], 0.0);
  // fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[0][5], 0.0);
  // fprintf(file_inv_time, "%lf %lf\n", simulation[0][6], 0.0);

	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[i][3], sim_time);
  //     fprintf(file_obj_time, "%lf %lf\n", simulation[i][4], sim_time);
  //     fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[i][5], sim_time);
  //     fprintf(file_inv_time, "%lf %lf\n", simulation[i][6], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }



  // fclose(file_error);
  // fclose(file_time);
  // fclose(file_abs_error);
  // fclose(file_abs_time);
  // fclose(file_abs_normalized);
  // fclose(file_abs_normalized_time);
  // fclose(file_obj);
  // fclose(file_obj_time);
  // fclose(file_obj_normalized);
  // fclose(file_obj_normalized_time);
  // fclose(file_inv);
  // fclose(file_inv_time);



  // initialize(opt);
  // total_time = 0;
  // algo_name = "FW_syn";
  // cout << "start " +algo_name  + " algorithm..." << endl;
  // //multiplicative error and time--------------------------
  // file_name_cpp = folder_name + algo_name + mul_name;
  // file_name_c = file_name_cpp.c_str();
  // file_error = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + mul_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_time = fopen(file_name_c, "w");
  // //absolute error and time--------------------------
  // file_name_cpp = folder_name + algo_name + abs_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_error = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + abs_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_time = fopen(file_name_c, "w");
  // //normalized absolute error and time--------------------------
  // file_name_cpp = folder_name + algo_name + abs_normalized_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_normalized = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_normalized_time = fopen(file_name_c, "w");
  // //number of inversions and time-------------------------
  // file_name_cpp = folder_name + algo_name + inv_name;
  // file_name_c = file_name_cpp.c_str();
  // file_inv = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + inv_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_inv_time = fopen(file_name_c, "w");
  // //objective function and time/ normalized objtive function and time.-------------------------
  // file_name_cpp = folder_name + algo_name + obj_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_norm_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_normalized = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_time = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_norm_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_normalized_time = fopen(file_name_c, "w");
  
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;



  // error = multiplicative(opt);
  // abs_error = absolute(opt);
  // obj = obj_function(opt);

  // //printf("the error is %lf\n", error);
  // simulation[0][0] = error; //multiplicative error
  // simulation[0][1] = (double)(0);
  // simulation[0][2] = abs_error; //absolute error
  // simulation[0][3] = abs_error / real_norm;
  // simulation[0][4] = obj;
  // simulation[0][5] = obj - real_obj;
  // fprintf(file_error, "%d %lf\n", (0), error);
  // fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
  // fprintf(file_abs_error, "%d %lf\n", (0), simulation[0][2]);
  // for (int i1 =0; i1 < opt->n; i1++){
  //   Mark_nodes[i1] = opt->hypernodes[i1];
  // }
  // qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  // inv = num_inversions(Mark_nodes, 0, opt->n-1);
  // simulation[0][6] = inv;
  // fprintf(file_inv, "%d %d\n", (0), inv);
  // fprintf(file_obj_normalized, "%d %lf\n", (0), simulation[0][5]);
  // fprintf(file_abs_normalized, "%d %lf\n", (0), simulation[0][3]);
  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
  //   //onepass_PR_plus_normal(opt, last_x, z, &lambda);
	// 	//onepass_PR_plus_new_2_normal_change(opt, last_x, z, last_gradient);
  //   onepass_FW_syn(opt);
  //   //onepass_PR_plus_new_2(opt, last_x, z);
  //   time2 = getCurrentTimeInMilliseconds();
  //     //test
  //   // for (int i1 = 0; i1 < opt->n; i1++){
  //   //   printf("r = %lf\n", opt->hypernodes[i1].r);
  //   // }
  //   // scanf("%d", &rep);
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative(opt);
  //     abs_error = absolute(opt);
  //     obj = obj_function(opt);

  //     //printf("the error is %lf\n", error);
  //     simulation[i][0] = error; //multiplicative error
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = abs_error; //absolute error
  //     simulation[i][3] = abs_error / real_norm;
  //     simulation[i][4] = obj;
  //     simulation[i][5] = obj - real_obj;
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj);
  //     fprintf(file_abs_error, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     inv = num_inversions(Mark_nodes, 0, opt->n-1);
  //     simulation[i][6] = inv;
  //     fprintf(file_inv, "%d %d\n", (i+1), inv);
  //     fprintf(file_obj_normalized, "%d %lf\n", (i+1), simulation[i][5]);
  //     fprintf(file_abs_normalized, "%d %lf\n", (i+1), simulation[i][3]);
  //   }
	// }


  // average_time = (double)total_time / (double)rep;
  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_abs_time, "%lf %lf\n", simulation[0][2], 0.0);
  // fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[0][3], 0.0);
  // fprintf(file_obj_time, "%lf %lf\n", simulation[0][4], 0.0);
  // fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[0][5], 0.0);
  // fprintf(file_inv_time, "%lf %lf\n", simulation[0][6], 0.0);

	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[i][3], sim_time);
  //     fprintf(file_obj_time, "%lf %lf\n", simulation[i][4], sim_time);
  //     fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[i][5], sim_time);
  //     fprintf(file_inv_time, "%lf %lf\n", simulation[i][6], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }



  // fclose(file_error);
  // fclose(file_time);
  // fclose(file_abs_error);
  // fclose(file_abs_time);
  // fclose(file_abs_normalized);
  // fclose(file_abs_normalized_time);
  // fclose(file_obj);
  // fclose(file_obj_time);
  // fclose(file_obj_normalized);
  // fclose(file_obj_normalized_time);
  // fclose(file_inv);
  // fclose(file_inv_time);


  // initialize(opt);
	// for (int i1 = 0; i1 < opt->e; i1++){
	// 	last_x[i1].a = opt->edges[i1].a;
  //   z[i1].a = opt->edges[i1].a;
  //   y[i1].a = opt->edges[i1].a;
	// 	//last_x[i].b = 1 - opt->edges[i].a;
	// }
  // total_time = 0;
  // algo_name = "FISTA";
  // cout << "start " +algo_name  + " algorithm..." << endl;
  // //multiplicative error and time--------------------------
  // file_name_cpp = folder_name + algo_name + mul_name;
  // file_name_c = file_name_cpp.c_str();
  // file_error = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + mul_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_time = fopen(file_name_c, "w");
  // //absolute error and time--------------------------
  // file_name_cpp = folder_name + algo_name + abs_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_error = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + abs_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_time = fopen(file_name_c, "w");
  // //normalized absolute error and time--------------------------
  // file_name_cpp = folder_name + algo_name + abs_normalized_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_normalized = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + abs_normalized_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_abs_normalized_time = fopen(file_name_c, "w");
  // //number of inversions and time-------------------------
  // file_name_cpp = folder_name + algo_name + inv_name;
  // file_name_c = file_name_cpp.c_str();
  // file_inv = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + inv_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_inv_time = fopen(file_name_c, "w");
  // //objective function and time/ normalized objtive function and time.-------------------------
  // file_name_cpp = folder_name + algo_name + obj_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_norm_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_normalized = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_time = fopen(file_name_c, "w");
  // file_name_cpp = folder_name + algo_name + obj_norm_time_name;
  // file_name_c = file_name_cpp.c_str();
  // file_obj_normalized_time = fopen(file_name_c, "w");
  
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;



  // error = multiplicative(opt);
  // abs_error = absolute(opt);
  // obj = obj_function(opt);

  // //printf("the error is %lf\n", error);
  // simulation[0][0] = error; //multiplicative error
  // simulation[0][1] = (double)(0);
  // simulation[0][2] = abs_error; //absolute error
  // simulation[0][3] = abs_error / real_norm;
  // simulation[0][4] = obj;
  // simulation[0][5] = obj - real_obj;
  // fprintf(file_error, "%d %lf\n", (0), error);
  // fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
  // fprintf(file_abs_error, "%d %lf\n", (0), simulation[0][2]);
  // for (int i1 =0; i1 < opt->n; i1++){
  //   Mark_nodes[i1] = opt->hypernodes[i1];
  // }
  // qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  // inv = num_inversions(Mark_nodes, 0, opt->n-1);
  // simulation[0][6] = inv;
  // fprintf(file_inv, "%d %d\n", (0), inv);
  // fprintf(file_obj_normalized, "%d %lf\n", (0), simulation[0][5]);
  // fprintf(file_abs_normalized, "%d %lf\n", (0), simulation[0][3]);
  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
  //   //onepass_PR_plus_normal(opt, last_x, z, &lambda);
	// 	//onepass_PR_plus_new_2_normal_change(opt, last_x, z, last_gradient);
  //   onepass_FISTA_normal(opt, max_degree, z, y, last_x);
  //   //onepass_PR_plus_new_2(opt, last_x, z);
  //   time2 = getCurrentTimeInMilliseconds();
  //     //test
  //   // for (int i1 = 0; i1 < opt->n; i1++){
  //   //   printf("r = %lf\n", opt->hypernodes[i1].r);
  //   // }
  //   // scanf("%d", &rep);
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative(opt);
  //     abs_error = absolute(opt);
  //     obj = obj_function(opt);

  //     //printf("the error is %lf\n", error);
  //     simulation[i][0] = error; //multiplicative error
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = abs_error; //absolute error
  //     simulation[i][3] = abs_error / real_norm;
  //     simulation[i][4] = obj;
  //     simulation[i][5] = obj - real_obj;
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj);
  //     fprintf(file_abs_error, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     inv = num_inversions(Mark_nodes, 0, opt->n-1);
  //     simulation[i][6] = inv;
  //     fprintf(file_inv, "%d %d\n", (i+1), inv);
  //     fprintf(file_obj_normalized, "%d %lf\n", (i+1), simulation[i][5]);
  //     fprintf(file_abs_normalized, "%d %lf\n", (i+1), simulation[i][3]);
  //   }
	// }


  // average_time = (double)total_time / (double)rep;
  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_abs_time, "%lf %lf\n", simulation[0][2], 0.0);
  // fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[0][3], 0.0);
  // fprintf(file_obj_time, "%lf %lf\n", simulation[0][4], 0.0);
  // fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[0][5], 0.0);
  // fprintf(file_inv_time, "%lf %lf\n", simulation[0][6], 0.0);

	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     fprintf(file_abs_normalized_time, "%lf %lf\n", simulation[i][3], sim_time);
  //     fprintf(file_obj_time, "%lf %lf\n", simulation[i][4], sim_time);
  //     fprintf(file_obj_normalized_time, "%lf %lf\n", simulation[i][5], sim_time);
  //     fprintf(file_inv_time, "%lf %lf\n", simulation[i][6], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }



  // fclose(file_error);
  // fclose(file_time);
  // fclose(file_abs_error);
  // fclose(file_abs_time);
  // fclose(file_abs_normalized);
  // fclose(file_abs_normalized_time);
  // fclose(file_obj);
  // fclose(file_obj_time);
  // fclose(file_obj_normalized);
  // fclose(file_obj_normalized_time);
  // fclose(file_inv);
  // fclose(file_inv_time);





  // return 1;

//-----------------------------------------------------------------------------------------


  // initialize(opt);
  


  // printf("start Greedy++ algorithm...\n");
  // file_error = fopen("data1/Greedy.txt", "w");
  // file_Abserror = fopen("data1/Greedy_abs.txt", "w");
  // file_inv = fopen("data1/Greedy_inv.txt", "w");
  // file_time = fopen("data1/Greedy_time.txt", "w");
  // file_Abs_time = fopen("data1/Greedy_Abs_time.txt", "w");
  // file_obj = fopen("data1/Greedy_obj.txt", "w");
  // total_time = 0;
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;
  

  //   for (int i1 =0; i1 < opt->n; i1++){
  //     Mark_nodes[i1].r = opt->hypernodes[i1].r;
  //   }

  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);

  //   time1 = getCurrentTimeInMilliseconds();
	// 	onepass_Greedypp(opt);
  //   time2 = getCurrentTimeInMilliseconds();
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
     
  //   total_time = total_time + time2 - time1;
  //   //printf("iteration %d finished\n", i+1);
    
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative_greedy(opt, (i+1));
  //     simulation[i][0] = error;
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = absolute_greedy(opt, i+1);
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj_function2(opt, i+1));
  //     fprintf(file_Abserror, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     num = num_inversions(Mark_nodes, 0, opt->n-1);
  //     fprintf(file_inv, "%d %d\n", (i+1), num);
  //   }
	// }
  // average_time = (double)total_time / (double)rep;

	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }

	// }
  // fclose(file_Abserror);
  // fclose(file_error);
  // fclose(file_inv);
  // fclose(file_time);
  // fclose(file_Abs_time);
  // fclose(file_obj);
  

  // printf("start Elist++ algorithm...\n");
  // initialize(opt);
  // if (opt->iter == 0){
  //   for (int i = 0; i < opt->n;i++){
  //     opt->hypernodes[i].r = 0;
  //   }
  // }
  // file_error = fopen("data1/Elist.txt", "w");
  // file_Abserror = fopen("data1/Elist_abs.txt", "w");
  // file_inv = fopen("data1/Elist_inv.txt", "w");
  // file_time = fopen("data1/Elist_time.txt", "w");
  // file_Abs_time = fopen("data1/Elist_Abs_time.txt", "w");
  // file_obj = fopen("data1/Elist_obj.txt", "w");
  // total_time = 0;
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;
  

  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
	// 	onepass_Elistpp(opt);
  //   time2 = getCurrentTimeInMilliseconds();
     
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
    
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative2(opt, (i+1));
  //     simulation[i][0] = error;
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = absolute2(opt, i+1);
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj_function2(opt, i+1));
  //     fprintf(file_Abserror, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     num = num_inversions(Mark_nodes, 0, opt->n-1);
  //     fprintf(file_inv, "%d %d\n", (i+1), num);
  //   }
	// }
  // average_time = (double)total_time / (double)rep;

	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }

  // fclose(file_Abserror);
  // fclose(file_error);
  // fclose(file_inv);
  // fclose(file_time);
  // fclose(file_Abs_time);
  // fclose(file_obj);
  // initialize(opt);


  // printf("start PR algorithm...\n");
  // file_error = fopen("data1/PR.txt", "w");
  // file_Abserror = fopen("data1/PR_abs.txt", "w");
  // file_inv = fopen("data1/PR_inv.txt", "w");
  // file_time = fopen("data1/PR_time.txt", "w");
  // file_Abs_time = fopen("data1/PR_Abs_time.txt", "w");
  // file_obj = fopen("data1/PR_obj.txt", "w");
  // total_time = 0;
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;


  // error = multiplicative(opt);
  // //printf("the error is %lf\n", error);
  // simulation[0][0] = error;
  // simulation[0][1] = (double)(0);
  // simulation[0][2] = absolute(opt);
  // fprintf(file_error, "%d %lf\n", (0), error);
  // fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
  // fprintf(file_Abserror, "%d %lf\n", (0), simulation[0][2]);
  // for (int i1 =0; i1 < opt->n; i1++){
  //   Mark_nodes[i1] = opt->hypernodes[i1];
  // }
  // qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  // num = num_inversions(Mark_nodes, 0, opt->n-1);
  // fprintf(file_inv, "%d %d\n", (0), num);


  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
	// 	onepass_PR(opt);
  //   time2 = getCurrentTimeInMilliseconds();
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
    
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative(opt);
  //     simulation[i][0] = error;
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = absolute(opt);
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj_function(opt));
  //     fprintf(file_Abserror, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     num = num_inversions(Mark_nodes, 0, opt->n-1);
  //     fprintf(file_inv, "%d %d\n", (i+1), num);
  //   }
	// }
  // average_time = (double)total_time / (double)rep;
  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_Abs_time, "%lf %lf\n", simulation[0][2], 0.0);
	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }
  // fclose(file_error);
  // fclose(file_inv);
  // fclose(file_time);
  // fclose(file_Abserror);
  // fclose(file_Abs_time);
  // fclose(file_obj);
  // initialize(opt);
  // total_time = 0;

  // printf("start PR_plus algorithm...\n");
  // file_error = fopen("data1/PR_plus_lin.txt", "w");
  // file_Abserror = fopen("data1/PR_plus_abs_lin.txt", "w");
  // file_inv = fopen("data1/PR_plus_inv_lin.txt", "w");
  // file_time = fopen("data1/PR_plus_time_lin.txt", "w");
  // file_Abs_time = fopen("data1/PR_plus_Abs_time_lin.txt", "w");
  // file_obj = fopen("data1/PR_plus_obj_lin.txt", "w");
  
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;


	// for (int i1 = 0; i1 < opt->e; i1++){
	// 	last_x[i1].a = opt->edges[i1].a;
  //   z[i1].a = opt->edges[i1].a;
	// 	//last_x[i].b = 1 - opt->edges[i].a;
	// }

  // error = multiplicative(opt);
  // //printf("the error is %lf\n", error);
  // simulation[0][0] = error;
  // simulation[0][1] = (double)(0);
  // simulation[0][2] = absolute(opt);
  // fprintf(file_error, "%d %lf\n", (0), error);
  // fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
  // fprintf(file_Abserror, "%d %lf\n", (0), simulation[0][2]);
  // for (int i1 =0; i1 < opt->n; i1++){
  //   Mark_nodes[i1] = opt->hypernodes[i1];
  // }
  // qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  // num = num_inversions(Mark_nodes, 0, opt->n-1);
  // fprintf(file_inv, "%d %d\n", (0), num);

  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
  //   //onepass_PR_plus_normal(opt, last_x, z, &lambda);
	// 	onepass_PR_plus_new_2_normal(opt, last_x, z);
  //   //onepass_PR_plus_product_normal(opt, last_x, z);
  //   //onepass_PR_plus_new_2(opt, last_x, z);
  //   time2 = getCurrentTimeInMilliseconds();
  //     //test
  //   // for (int i1 = 0; i1 < opt->n; i1++){
  //   //   printf("r = %lf\n", opt->hypernodes[i1].r);
  //   // }
  //   // scanf("%d", &rep);
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative(opt);
  //     simulation[i][0] = error;
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = absolute(opt);
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj_function(opt));
  //     fprintf(file_Abserror, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     num = num_inversions(Mark_nodes, 0, opt->n-1);
  //     fprintf(file_inv, "%d %d\n", (i+1), num);
  //   }
	// }

  // average_time = (double)total_time / (double)rep;

  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_Abs_time, "%lf %lf\n", simulation[0][2], 0.0);
	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }



  // fclose(file_error);
  // fclose(file_inv);
  // fclose(file_time);
  // fclose(file_Abserror);
  // fclose(file_Abs_time);
  // fclose(file_obj);
  // initialize(opt);


  // total_time = 0;

  // printf("start PR_plus algorithm...\n");
  // file_error = fopen("data1/PR_plus_exp_pro.txt", "w");
  // file_Abserror = fopen("data1/PR_plus_abs_exp_pro.txt", "w");
  // file_inv = fopen("data1/PR_plus_inv_exp_pro.txt", "w");
  // file_time = fopen("data1/PR_plus_time_exp_pro.txt", "w");
  // file_Abs_time = fopen("data1/PR_plus_Abs_time_exp_pro.txt", "w");
  // file_obj = fopen("data1/PR_plus_obj_exp_pro.txt", "w");
  
  // // double error = 0;
  // // struct timeval Time1;
  // // struct timeval Time2;


	// for (int i1 = 0; i1 < opt->e; i1++){
	// 	last_x[i1].a = opt->edges[i1].a;
  //   z[i1].a = opt->edges[i1].a;
	// 	//last_x[i].b = 1 - opt->edges[i].a;
	// }

  // error = multiplicative(opt);
  // //printf("the error is %lf\n", error);
  // simulation[0][0] = error;
  // simulation[0][1] = (double)(0);
  // simulation[0][2] = absolute(opt);
  // fprintf(file_error, "%d %lf\n", (0), error);
  // fprintf(file_obj, "%d %lf\n", (0), obj_function(opt));
  // fprintf(file_Abserror, "%d %lf\n", (0), simulation[0][2]);
  // for (int i1 =0; i1 < opt->n; i1++){
  //   Mark_nodes[i1] = opt->hypernodes[i1];
  // }
  // qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  // num = num_inversions(Mark_nodes, 0, opt->n-1);
  // fprintf(file_inv, "%d %d\n", (0), num);

  // //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	// for (i=0;i<rep;i++){
	// 	//printf("%u\n",i);
	// 	//onepass(opt);
  //   //onepass_PR_plus(opt, last_x, z, &lambda);
	// 	//onepass_PR(opt);
  //   //onepass_Greedypp(opt);
  //   //onepass_Elistpp(opt);
  //   time1 = getCurrentTimeInMilliseconds();
  //   //onepass_PR_plus_normal(opt, last_x, z, &lambda);
	// 	//onepass_PR_plus_new_2_normal(opt, last_x, z);
  //   onepass_PR_plus_product_normal_2(opt, last_x, z);
  //   //onepass_PR_plus_new_2(opt, last_x, z);
  //   time2 = getCurrentTimeInMilliseconds();
  //     //test
  //   // for (int i1 = 0; i1 < opt->n; i1++){
  //   //   printf("r = %lf\n", opt->hypernodes[i1].r);
  //   // }
  //   // scanf("%d", &rep);
  //   total_time = total_time + time2 - time1;
  //   //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
  //   //onepass_FISTA(opt, max_degree, z, y, last_x);
  //   //printf("iteration %d finished\n", i+1);
  //   //printf("the error is %lf\n", error);
  //   if ((i+1)%num_points == 0){
  //     error = multiplicative(opt);
  //     simulation[i][0] = error;
  //     simulation[i][1] = (double)(i+1);
  //     simulation[i][2] = absolute(opt);
  //     fprintf(file_error, "%d %lf\n", (i+1), error);
  //     fprintf(file_obj, "%d %lf\n", (i+1), obj_function(opt));
  //     fprintf(file_Abserror, "%d %lf\n", (i+1), simulation[i][2]);
  //     for (int i1 =0; i1 < opt->n; i1++){
  //       Mark_nodes[i1] = opt->hypernodes[i1];
  //     }
  //     qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  //     num = num_inversions(Mark_nodes, 0, opt->n-1);
  //     fprintf(file_inv, "%d %d\n", (i+1), num);
  //   }
	// }

  // average_time = (double)total_time / (double)rep;

  // fprintf(file_time, "%lf %lf\n", simulation[0][0], 0.0);
  // fprintf(file_Abs_time, "%lf %lf\n", simulation[0][2], 0.0);
	// for (i=0;i<rep;i++){
  //   if ((i+1)%num_points == 0){
  //     double sim_time = simulation[i][1]*average_time;
  //     fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
  //     fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
  //     if (sim_time < 0){
  //       printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
  //     }
  //   }
	// }



  // fclose(file_error);
  // fclose(file_inv);
  // fclose(file_time);
  // fclose(file_Abserror);
  // fclose(file_Abs_time);
  // fclose(file_obj);
  // initialize(opt);


  //----------------------------------------------------------------




/*
//------------------------------------------------------------------------------------------------------------

  file_error = fopen("data1/PR_plus_1.txt", "w");
  file_Abserror = fopen("data1/PR_plus_abs_1.txt", "w");
  file_inv = fopen("data1/PR_plus_inv_1.txt", "w");
  file_time = fopen("data1/PR_plus_time_1.txt", "w");
  file_Abs_time = fopen("data1/PR_plus_Abs_time_1.txt", "w");
  
  // double error = 0;
  // struct timeval Time1;
  // struct timeval Time2;

  error = multiplicative(opt);
  //printf("the error is %lf\n", error);
  simulation[0][0] = error;
  simulation[0][1] = (double)(0);
  simulation[0][2] = absolute(opt);
  fprintf(file_error, "%d %lf\n", (0), error);
  fprintf(file_Abserror, "%d %lf\n", (0), simulation[0][2]);
  for (int i1 =0; i1 < opt->n; i1++){
    Mark_nodes[i1] = opt->hypernodes[i1];
  }
  qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  num = num_inversions(Mark_nodes, 0, opt->n-1);
  fprintf(file_inv, "%d %d\n", (0), num);


	for (int i1 = 0; i1 < opt->e; i1++){
		last_x[i1].a = opt->edges[i1].a;
    z[i1].a = opt->edges[i1].a;
		//last_x[i].b = 1 - opt->edges[i].a;
	}
  

  //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, z, &lambda);
		//onepass_PR(opt);
    //onepass_Greedypp(opt);
    //onepass_Elistpp(opt);
    time1 = getCurrentTimeInMilliseconds();
		onepass_PR_plus_new_1(opt, last_x, z);
    time2 = getCurrentTimeInMilliseconds();
     
    total_time = total_time + time2 - time1;
    //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
    //printf("iteration %d finished\n", i+1);
    update_density(opt);
    error = multiplicative(opt);
    //printf("the error is %lf\n", error);
    if ((i+1)%10 == 0){
      simulation[i][0] = error;
      simulation[i][1] = (double)(i+1);
      fprintf(file_error, "%d %lf\n", (i+1), error);
      fprintf(file_Abserror, "%d %lf\n", (i+1), absolute(opt));
      for (int i1 =0; i1 < opt->n; i1++){
        Mark_nodes[i1] = opt->hypernodes[i1];
      }
      qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
      int num = num_inversions(Mark_nodes, 0, opt->n-1);
      fprintf(file_inv, "%d %d\n", (i+1), num);
    }
	}

  average_time = (double)total_time / (double)rep;

	for (i=0;i<rep;i++){
    if ((i+1)%10 == 0){
      double sim_time = simulation[i][1]*average_time;
      fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
      fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
      if (sim_time < 0){
        printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
      }
      //test
      // for (int i1 = 0; i1 < opt->n; i1++){
      //   printf("r = %lf\n", opt->hypernodes[i1].r);
      // }
      // scanf("%d", &rep);
    }
	}




  fclose(file_error);
  fclose(file_inv);
  fclose(file_time);
  fclose(file_Abserror);
  initialize(opt);



  file_error = fopen("data1/PR_plus_2.txt", "w");
  file_Abserror = fopen("data1/PR_plus_abs_2.txt", "w");
  file_inv = fopen("data1/PR_plus_inv_2.txt", "w");
  file_time = fopen("data1/PR_plus_time_2.txt", "w");
  file_Abs_time = fopen("data1/PR_plus_Abs_time_2.txt", "w");
  
  // double error = 0;
  // struct timeval Time1;
  // struct timeval Time2;

	for (int i1 = 0; i1 < opt->e; i1++){
		last_x[i1].a = opt->edges[i1].a;
    z[i1].a = opt->edges[i1].a;
		//last_x[i].b = 1 - opt->edges[i].a;
	}

  error = multiplicative(opt);
  //printf("the error is %lf\n", error);
  simulation[0][0] = error;
  simulation[0][1] = (double)(0);
  simulation[0][2] = absolute(opt);
  fprintf(file_error, "%d %lf\n", (0), error);
  fprintf(file_Abserror, "%d %lf\n", (0), simulation[0][2]);
  for (int i1 =0; i1 < opt->n; i1++){
    Mark_nodes[i1] = opt->hypernodes[i1];
  }
  qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
  num = num_inversions(Mark_nodes, 0, opt->n-1);
  fprintf(file_inv, "%d %d\n", (0), num);



  //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, z, &lambda);
		//onepass_PR(opt);
    //onepass_Greedypp(opt);
    //onepass_Elistpp(opt);
    time1 = getCurrentTimeInMilliseconds();
		onepass_PR_plus_new_2_normal(opt, last_x, z);
    time2 = getCurrentTimeInMilliseconds();
     
    total_time = total_time + time2 - time1;
    //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
    //printf("iteration %d finished\n", i+1);
    update_density(opt);
    error = multiplicative(opt);
    //printf("the error is %lf\n", error);
    if ((i+1)%10 == 0){
      simulation[i][0] = error;
      simulation[i][1] = (double)(i+1);
      fprintf(file_error, "%d %lf\n", (i+1), error);
      fprintf(file_Abserror, "%d %lf\n", (i+1), absolute(opt));
      for (int i1 =0; i1 < opt->n; i1++){
        Mark_nodes[i1] = opt->hypernodes[i1];
      }
      qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
      int num = num_inversions(Mark_nodes, 0, opt->n-1);
      fprintf(file_inv, "%d %d\n", (i+1), num);
    }
	}

  average_time = (double)total_time / (double)rep;

	for (i=0;i<rep;i++){
    if ((i+1)%10 == 0){
      double sim_time = simulation[i][1]*average_time;
      fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
      fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
      if (sim_time < 0){
        printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
      }
      //test
      // for (int i1 = 0; i1 < opt->n; i1++){
      //   printf("r = %lf\n", opt->hypernodes[i1].r);
      // }
      // scanf("%d", &rep);
    }
	}




  fclose(file_error);
  fclose(file_inv);
  fclose(file_time);
  fclose(file_Abserror);
  initialize(opt);



  file_error = fopen("data1/PR_plus_1.txt", "w");
  file_Abserror = fopen("data1/PR_plus_abs_1.txt", "w");
  file_inv = fopen("data1/PR_plus_inv_1.txt", "w");
  file_time = fopen("data1/PR_plus_time_1.txt", "w");
  file_Abs_time = fopen("data1/PR_plus_Abs_time_1.txt", "w");
  
  // double error = 0;
  // struct timeval Time1;
  // struct timeval Time2;
  

  //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, z, &lambda);
		//onepass_PR(opt);
    //onepass_Greedypp(opt);
    //onepass_Elistpp(opt);
    time1 = getCurrentTimeInMilliseconds();
		onepass_PR_plus_new_1(opt, last_x, z);
    time2 = getCurrentTimeInMilliseconds();
     
    total_time = total_time + time2 - time1;
    //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
    //printf("iteration %d finished\n", i+1);
    update_density(opt);
    error = multiplicative(opt);
    //printf("the error is %lf\n", error);
    if ((i+1)%10 == 0){
      simulation[i][0] = error;
      simulation[i][1] = (double)(i+1);
      fprintf(file_error, "%d %lf\n", (i+1), error);
      fprintf(file_Abserror, "%d %lf\n", (i+1), absolute(opt));
      for (int i1 =0; i1 < opt->n; i1++){
        Mark_nodes[i1] = opt->hypernodes[i1];
      }
      qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
      int num = num_inversions(Mark_nodes, 0, opt->n-1);
      fprintf(file_inv, "%d %d\n", (i+1), num);
    }
	}

  average_time = (double)total_time / (double)rep;

	for (i=0;i<rep;i++){
    if ((i+1)%10 == 0){
      double sim_time = simulation[i][1]*average_time;
      fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
      fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
      if (sim_time < 0){
        printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
      }
      //test
      // for (int i1 = 0; i1 < opt->n; i1++){
      //   printf("r = %lf\n", opt->hypernodes[i1].r);
      // }
      // scanf("%d", &rep);
    }
	}




  fclose(file_error);
  fclose(file_inv);
  fclose(file_time);
  fclose(file_Abserror);
  initialize(opt);



  file_error = fopen("data1/PR_plus_2.txt", "w");
  file_Abserror = fopen("data1/PR_plus_abs_2.txt", "w");
  file_inv = fopen("data1/PR_plus_inv_2.txt", "w");
  file_time = fopen("data1/PR_plus_time_2.txt", "w");
  file_Abs_time = fopen("data1/PR_plus_Abs_time_2.txt", "w");
  
  // double error = 0;
  // struct timeval Time1;
  // struct timeval Time2;
  

  //printf("Current time: %ld seconds and %ld microseconds\n", currentTime.tv_sec, currentTime.tv_usec);
	for (i=0;i<rep;i++){
		//printf("%u\n",i);
		//onepass(opt);
    //onepass_PR_plus(opt, last_x, z, &lambda);
		//onepass_PR(opt);
    //onepass_Greedypp(opt);
    //onepass_Elistpp(opt);
    time1 = getCurrentTimeInMilliseconds();
		onepass_PR_plus_new_2(opt, last_x, z);
    time2 = getCurrentTimeInMilliseconds();
     
    total_time = total_time + time2 - time1;
    //fprintf(file_time, "%d %ld\n", (i+1), Time2.tv_usec -  Time1.tv_usec);
    //onepass_FISTA(opt, max_degree, z, y, last_x);
    //printf("iteration %d finished\n", i+1);
    update_density(opt);
    error = multiplicative(opt);
    //printf("the error is %lf\n", error);
    if ((i+1)%10 == 0){
      simulation[i][0] = error;
      simulation[i][1] = (double)(i+1);
      fprintf(file_error, "%d %lf\n", (i+1), error);
      fprintf(file_Abserror, "%d %lf\n", (i+1), absolute(opt));
      for (int i1 =0; i1 < opt->n; i1++){
        Mark_nodes[i1] = opt->hypernodes[i1];
      }
      qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);
      int num = num_inversions(Mark_nodes, 0, opt->n-1);
      fprintf(file_inv, "%d %d\n", (i+1), num);
    }
	}

  average_time = (double)total_time / (double)rep;

	for (i=0;i<rep;i++){
    if ((i+1)%10 == 0){
      double sim_time = simulation[i][1]*average_time;
      fprintf(file_time, "%lf %lf\n", simulation[i][0], sim_time);
      fprintf(file_Abs_time, "%lf %lf\n", simulation[i][2], sim_time);
      if (sim_time < 0){
        printf("total time = %ld, rep = %d, average_time = %lf, sim_time = %lf\n", total_time, rep, average_time, sim_time);
      }
      //test
      for (int i1 = 0; i1 < opt->n; i1++){
        printf("r = %lf\n", opt->hypernodes[i1].r);
      }
      scanf("%d", &rep);
    }
	}




  fclose(file_error);
  fclose(file_inv);
  fclose(file_time);
  fclose(file_Abserror);
  initialize(opt);










  if (flag_pp == 1){
    for (int i=0; i < opt->n; i++){
      opt->hypernodes[i].r = opt->hypernodes[i].r / rep;
    }
  }



  //node *Mark_nodes = (node*)malloc(opt->n * sizeof(node));
  for (int i =0; i < opt->n; i++){
    Mark_nodes[i] = opt->hypernodes[i];
  }

  qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes2);

  
  int num = num_inversions(Mark_nodes, 0, opt->n-1);
  printf("laksdjflkd\n");

  printf("the number of inversions is = %d\n", num);
  update_density(opt);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
  // Out_put_density(opt, "output_PR1.txt");
  // Out_put_weights(opt, "output_PR1w.txt");
  printf("%d, %d\n", opt->m, opt->n);
  node *Mark_nodes = (node*)malloc(opt->n * sizeof(node));
  for (int i =0; i < opt->n; i++){
    Mark_nodes[i] = opt->hypernodes[i];
  }
  //qsort(Mark_nodes,opt->n,sizeof(node),compare_nodes);

  t3=time(NULL);
  int res = verify(opt, Mark_nodes);
  printf("verify, %d", res);
  t4 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n",(t4-t3)/3600,((t4-t3)%3600)/60,((t4-t3)%60));

  FILE *exact = fopen("Exact.txt", "w");

  for (int i = 0; i < opt->n; i++){
    fprintf(exact, "%d %.20f\n", i, opt->hypernodes[i].r_correct);
    if (opt->hypernodes[i].old_id != i) printf("strange");
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

*/

