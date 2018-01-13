#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <set>

/*
  gen.cpp 
  
  Usage: 
  gen <number of point> <number of edge>
  
  output format:
  n m n1 n2
  <node_weight> <node_type>
  ...
  <edge_from> <edge_to>
  ...
  <I_to_O_weight>
  
  n = |V|  size of the node set V
  m = |E|  size of the edge set E
  n1 = |I| size of the set I
  n2 = |O| size of the set O
  
 */

#define MAX_WEIGHT 10000
using namespace std;


struct node{
  int p;
  node *next;
  node(int p = 0, node *next = NULL): p(p), next(next) {} 
};

int get_rand_weight(){
  return rand() % MAX_WEIGHT;
}

int main(int argc, char *argv[]){
  if(argc != 5){
    printf("Usage: gen <number of point> <number of edge> <number of point in I> <number of point in O>\n");
    exit(0);
  }

  //freopen("data.txt","w",stdout);
  srand(time(NULL));
  
  int n = atoi(argv[1]);
  int m = atoi(argv[2]);
  int n1 = atoi(argv[3]);
  int n2 = atoi(argv[4]);
  
  int *node_weight = new int[n];
  node **head = new node*[n];

  static set<pair<int,int> > S;
  for(int i = 0; i < n; ++ i){
    node_weight[i] = get_rand_weight(); 
    head[i] = NULL;
  }

  int num = 0;
  for(int i = n1; i < n - n2; ++ i){
    int x = rand() % n1;
    while(S.find(make_pair(x, i)) != S.end()) x = rand() % n1 ;
    node *p = new node(i, head[x]);
    head[x] = p;
    S.insert(make_pair(x, i));
    ++ num;
    x = n - n2 + rand() % n2;
    while(S.find(make_pair(i, x)) != S.end()) x = n - n2 + rand() % n2 ;
    p = new node(x, head[i]);
    head[i] = p;
    S.insert(make_pair(i , x));
    ++ num;
  }

  for(int i = 0 ; i < n1; ++ i){
    int x = rand() % (n - n1) + n1;
    while(S.find(make_pair(i, x)) != S.end()) x = rand() % (n - n1) + n1;
    node *p = new node(x, head[i]);
    head[i] = p;
    S.insert(make_pair(i, x));
  }

  for(int i = 1 ; i <= n2; ++ i){
    int ii = n - i;
    int x = rand() % (n - n2);
    while(S.find(make_pair(x, ii)) != S.end()) x = rand() % (n - n2);
    node *p = new node(ii, head[x]);
    head[x] = p;
    S.insert(make_pair(x, ii));
  }
  
  for(int i = num; i < m ; ++ i){
    int x = rand() % (n - n2 - n1) + n1 , y = rand() % (n - n2 - n1) + n1; 
    while(x == y || S.find(make_pair(min(x,y), max(x, y))) != S.end())
      x = rand() % (n - n2 - n1) + n1 , y = rand() % (n - n2 - n1) + n1;
    if(x > y) swap(x, y);
    node *p = new node(y, head[x]);
    head[x] = p;
    S.insert(make_pair(x, y));
  }

  printf("%d %d %d %d\n",n,m,n1,n2);
  for(int i = 0; i < n; ++ i){
    int cmd = 0;
    if(i < n1) cmd = 1;
    else if(i >= n - n2) cmd = 2;
    printf("%d %d\n",node_weight[i],cmd);
  }

  for(int i = 0; i < n; ++ i)
    for(node *p = head[i]; p; p = p -> next)
      printf("%d %d\n",i, p -> p);
  //static set<pair<int,int>> Q;
  for(int i = 0; i < n1; ++ i, printf("\n"))
    for(int j = 0; j < n2; ++ j)
      printf("%d ",get_rand_weight());
  
  return 0;
}
