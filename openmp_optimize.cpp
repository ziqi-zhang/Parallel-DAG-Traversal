#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <vector>
#define pii pair<int,int>
#define MAXN 20000000
using namespace std;

static double get_time(){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec+(double)tv.tv_usec*1e-6;
}

struct node{
  int p;
  node *nxt;
  node(){
    p = 0;
    nxt = NULL;
  }
  //node(int p, node *next = NULL): p(p), next(next) {} 
};

struct node_info{
  int* top;
  node_info(){
    top = new int[10];
    memset(top, -1, sizeof(int) * 10);
  }
};

vector<int> I, O;
int main(int argc, char ** argv){
  int n, m, n1, n2;
  scanf("%d %d %d %d",&n,&m,&n1,&n2);

  int *node_weight = new int[n];
  int *node_type = new int[n];
  node **head = new node*[n];
  node **head_inv = new node*[n];
  
  for(int i = 0; i < n; ++ i){
    scanf("%d %d",node_weight + i, node_type + i);
    if(node_type[i] == 1)
      I.push_back(i);
    else
      if(node_type[i] == 2)
	O.push_back(i);
    head[i] = head_inv[i] = NULL;
  }

  node * nodebuf = new node[m * 2];

  
  for(int i = 0; i < m; ++ i){
    int u, v;
    scanf("%d %d", &u, &v);
    node *p = nodebuf + i;
    p -> p = v;
    p -> nxt = head[u];
    head[u] = p;

    p = nodebuf + m + i;
    p -> p = u;
    p -> nxt = head_inv[v];
    head_inv[v] = p;
  }

  int *weight_table = new int[n1 * n2];
  
  for(int i = 0; i < n1; ++ i)
    for(int j = 0; j < n2; ++ j)
      scanf("%d",&weight_table[i * n2 + j]);
  fclose(stdin);

  double t0 = get_time();

  omp_lock_t* lock = new omp_lock_t[n];
  omp_lock_t que_lock;

  
  for(int i = 0; i < n; ++ i)
    omp_init_lock(&lock[i]);
  omp_init_lock(&que_lock);
  
  int *du = new int[n];
  memset(du, 0, sizeof(int) * n);
  node_info * info = new node_info[n * n2];

  int nthreads = 20;
  
#pragma omp parallel for num_threads(nthreads)
  for(int i = 0; i < n; ++ i)
    for(node *p = head[i]; p; p = p -> nxt){
      omp_set_lock(&lock[p -> p]);
      du[p -> p] ++;
      omp_unset_lock(&lock[p -> p]);
    }

  int *que0 = new int[n];
  int *que1 = new int[n];
  int num0, num1;
  num0 = num1 = 0;
  for(int i = 0; i < n1; ++ i)
    que0[num0 ++] = I[i];
  
#pragma omp parallel for num_threads(nthreads)
  for(int i = 0; i < n1; ++ i)  
    for(int j = 0; j < n2; ++ j)
      info[I[i] * n2 + j].top[0] = weight_table[i * n2 + j];

  int round = 0;
  
  while(num0){
    num1 = 0;
    
    //int ii;
#pragma omp parallel for num_threads(nthreads)
    for(int ii = 0; ii < num0; ++ ii){
      int u = que0[ii];
      //cerr << "ok " << u << endl;
      for(node *i = head[u]; i; i = i -> nxt){
	int v = i -> p;
	omp_set_lock(&lock[v]);
	-- du[v];
	if(du[v] == 0){
	  omp_set_lock(&que_lock);
	  que1[num1 ++] = v;
	  //cerr << "ins " << v << endl; 
	  omp_unset_lock(&que_lock);
	}
	omp_unset_lock(&lock[v]);
      }
    }

    cerr << ++ round << ' ' << num0 << ' ' << num1 << endl;
#pragma omp parallel for num_threads(nthreads)
    for(int ii = 0; ii < num1; ++ ii){
      int v = que1[ii];
      
      for(node *p = head_inv[v]; p; p = p -> nxt){
	int u = p -> p;
	if(node_type[v] == 0){
	  int temp_top[10];	
	  for(int j = 0; j < n2; ++ j){
	    memset(temp_top, -1, sizeof(temp_top));
	    int h1 = 0, h2 = 0;
	    for(int k = 0; k < 10; ++ k){
	      if(info[v * n2 + j].top[h1] == -1 && info[u * n2 + j].top[h2] == -1) break;
	      if(info[v * n2 + j].top[h1] >= info[u * n2 + j].top[h2] + node_weight[v] || info[u * n2 + j].top[h2] == -1){
		temp_top[k] = info[v * n2 + j].top[h1];
		++ h1;
	      }else{
		temp_top[k] = info[u * n2 + j].top[h2] + node_weight[v];
		++ h2;
	      }
	    }
	    memcpy(info[v * n2 + j].top, temp_top, sizeof(temp_top));
	  }
	}else{
	  int temp_top[10];
	  for(int j = 0; j < n2; ++ j){
	    if(O[j] != v) continue;
	    memset(temp_top, -1, sizeof(temp_top));
	    int h1 = 0, h2 = 0;
	    for(int k = 0; k < 10; ++ k){
	      if(info[v * n2 + j].top[h1] == -1 && info[u * n2 + j].top[h2] == -1) break;
	      if(info[v * n2 + j].top[h1] >= info[u * n2 + j].top[h2]){
		temp_top[k] = info[v * n2 + j].top[h1];
		++ h1;
	      }else{
		temp_top[k] = info[u * n2 + j].top[h2];
		++ h2;
	      }
	    }
	    memcpy(info[v * n2 + j].top, temp_top, sizeof(temp_top));
	  }
	}
	
      }
      
    }
    swap(num0, num1);
    swap(que0, que1);
  }
  
  double t1 = get_time();
  cerr << "runtime:" << t1 - t0 << endl;
  freopen("openmp_result.txt","w",stdout);
      
  for(int i = 0; i < n2; ++ i){
    printf("%d ",O[i]);
    int x = O[i], cnt = 0;
    for(int j = 0; j < 10 && info[x * n2 + i].top[j] != -1; ++ j) ++ cnt;
    printf("%d ",cnt);
    for(int j = 0; j < cnt; ++ j)
      printf("%d ",info[x * n2 + i].top[j]);
    printf("\n");
  }
  
  return 0;
}

