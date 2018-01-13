#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <sys/time.h>
#include <time.h>
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
  
  for(int i = 0; i < n; ++ i){
    scanf("%d %d",node_weight + i, node_type + i);
    if(node_type[i] == 1)
      I.push_back(i);
    else
      if(node_type[i] == 2)
	O.push_back(i);
  }

  node * nodebuf = new node[m];
  
  for(int i = 0; i < m; ++ i){
    int u, v;
    scanf("%d %d", &u, &v);
    node *p = nodebuf + i;
    p -> p = v;
    p -> nxt = head[u];
    head[u] = p;
  }

  int *weight_table = new int[n1 * n2];
  
  for(int i = 0; i < n1; ++ i)
    for(int j = 0; j < n2; ++ j)
      scanf("%d",&weight_table[i * n2 + j]);
  fclose(stdin);

  double t0 = get_time();
  int *que = new int[n];
  int *du = new int[n];
  memset(du, 0, sizeof(int) * n);
  node_info * info = new node_info[n * n2];
  
  for(int i = 0; i < n; ++ i)
    for(node *p = head[i]; p; p = p -> nxt)
      du[p -> p] ++;

  int first, end;
  first = end = 0;
  for(int i = 0; i < n1; ++ i){
    que[end ++] = I[i];
    for(int j = 0; j < n2; ++ j)
      info[I[i] * n2 + j].top[0] = weight_table[i * n2 + j];
  }  
  /*
  for(int i = 0; i < n; ++ i, cerr << endl){
    for(int j = 0; j < 10; ++ j)
      cerr << info[i].top[j] << ' ';
  }
  return 0;*/
  while(first != end){
    if(first % 10 == 0) cerr << first << endl;
    int u = que[first++];
    for(node *i = head[u]; i; i = i -> nxt){
      int v = i -> p;
      -- du[v];
      if(du[v] == 0) que[end ++] = v;
      
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
  
  double t1 = get_time();
  cerr << "runtime:" << t1 - t0 << endl;
  freopen("serial_result.txt","w",stdout);
      
  for(int i = 0; i < n2; ++ i){
    printf("%d ",O[i]);
    int x = O[i], cnt = 0;
    for(int j = 0; j < 10 && info[x * n2 + i].top[j] != -1; ++ j) ++ cnt;
    printf("%d ",cnt);
    for(int j = 0; j < cnt; ++ j)
      printf("%d ",info[x * n2 + i].top[j]);
    printf("\n");
    // cerr << O[i] << endl;
    //for(int j = 0; j < cnt; ++ j) cerr << info[x].top[j] << ' '; cerr << endl;
  }
  
  return 0;
}

