#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <cuda_runtime.h>
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
};
//__device__ int my_block_count = 0;
//__device__ int my_block_count_1 = 0;

__global__ void expand(int num, int* inQue, int* edge_offset, int *edge_buf, int* du, int* lock){
  int tid = threadIdx.x, bid = blockIdx.x;
  for(int i = bid; i < num; i += gridDim.x){
    int u = inQue[i];
    int offset = edge_offset[u], offset_end = edge_offset[u + 1];
    for(int j = offset + tid; j < offset_end; j += blockDim.x){
      int y = edge_buf[j];
      //du[y] = bid;
      int d = atomicSub(&du[y], 1);
      if(d == 1)
	lock[y] = bid;
    }
  }
}

#define WARPSIZE 32
__global__ void count_inQue_offset(int num, int* inQue, int* edge_offset, int *edge_buf, int* du, int* lock, int* inQue_offset){
  int tid = threadIdx.x, bid = blockIdx.x;
  int accoffset = 0;
  for(int i = bid; i < num; i += gridDim.x){
    int u = inQue[i];
    int offset = edge_offset[u], offset_end = edge_offset[u + 1];
    int sum = 0;
    for(int j = offset + tid; j < offset_end; j += blockDim.x){
      int y = edge_buf[j];
      if(lock[y] == bid && du[y] == 0)
	++ sum;
    }
    
    __syncthreads();
    //reduce 
    for(int s = WARPSIZE / 2; s; s >>= 1)
      sum += __shfl_down(sum, s);
    
    if(tid == 0)
      accoffset += sum;

  }

  if(tid == 0)
    inQue_offset[bid] = accoffset;
}

__global__ void naive_scan(int num, int *data, int *sum){
  __shared__ int temp[WARPSIZE];
  
  int tid = threadIdx.x;
  int acc = 0;
  for(int i = 0; i < num; i += blockDim.x){
    int temp1;
    if(i + tid < num)
      temp[tid] = data[i + tid];
    else
      temp[tid] = 0;
    
    for(int d = 1; d < WARPSIZE; d <<= 1){
      if(tid >= d)
	temp1 = temp[tid - d];
      else
	temp1 = 0;
      __syncthreads();
      temp[tid] += temp1;
      __syncthreads();
    }

    if(i + tid < num)
      sum[i + tid] = temp[tid] + acc;
    __syncthreads();
    acc = sum[i + blockDim.x - 1];
  }
}


__global__ void enqueue(int num, int* inQue, int* outQue, int* edge_offset, int *edge_buf, int* lock, int *prefix_sum){
  __shared__ int acc ;
  
  int tid = threadIdx.x, bid = blockIdx.x;
  if(tid == 0){
    if(bid != 0)
      acc = prefix_sum[bid - 1];
    else
      acc = 0;
  }
  
  for(int i = bid; i < num; i += gridDim.x){
    int u = inQue[i];
    int offset = edge_offset[u], offset_end = edge_offset[u + 1];
    int sum = 0;
    for(int j = offset + tid; j < offset_end; j += blockDim.x){
      int y = edge_buf[j];
      if(lock[y] == bid)
	++ sum;
    }
    
    __syncthreads();
    //local scan
    int temp;
    for(int d = 1; d < WARPSIZE; d <<= 1){
      temp = __shfl_up(sum, d);
      if(tid >= d) sum += temp;
    }
    
    __syncthreads();
    int temp_acc = acc;
    __syncthreads();

    if(tid == 31){
      acc = sum;
    }

    sum = temp_acc + sum;
    
    for(int j = offset + tid; j < offset_end; j += blockDim.x){
      int y = edge_buf[j];
      if(lock[y] == bid)
	outQue[--sum] = y;
    }
    
    //__syncthreads();
    
  } 
}

#define nthreads 32
__global__ void update(int num, int n1, int n2,
		       int* node_weight,
		       int* node_type,
		       int* id,
		       int* outQue,
		       int* edge_inv_offset, int *edge_inv_buf,
		       int* top){
  __shared__ int shared_top[nthreads][10];
  
  int local_top[10];
  int temp_top[10];
  
  int tid = threadIdx.x, bid = blockIdx.x;
  for(int i = bid; i < num * n2; i += gridDim.x){
    for(int k = 0; k < 10; ++ k) local_top[k] = -1;
    int v = outQue[i / n2], o = i % n2;
    if(node_type[v] != 0 && id[o] != v) continue;
    int w = node_weight[v], offset = edge_inv_offset[v], offset_end = edge_inv_offset[v + 1];
    for(int j = offset + tid; j < offset_end; j += blockDim.x){
      
      int u = edge_inv_buf[j];
	  
      if(node_type[v] == 0){
	int h1 = 0, h2 = 0;
	int k;
	for(k = 0; k < 10; ++ k){
	  int top_u = top[(u * n2 + o) * 10 + h2];
	  if(top_u == -1 && local_top[h1] == -1) break;
	  if(top_u == -1 || local_top[h1] >= top_u + w)
	    temp_top[k] = local_top[h1], h1 ++ ;
	  else
	    temp_top[k] = top_u + w, h2 ++;
	}
	
	while(k < 10) temp_top[k++] = -1;
      }else{
	int h1 = 0, h2 = 0;
	int k;
	for(k = 0; k < 10; ++ k){
	  int top_u = top[(u * n2 + o) * 10 + h2];
	  if(top_u == -1 && local_top[h1] == -1) break;
	  if(top_u == -1 || local_top[h1] >= top_u)
	    temp_top[k] = local_top[h1], h1 ++ ;
	  else
	    temp_top[k] = top_u, h2 ++; 
	}
	
	while(k < 10) temp_top[k++] = -1;
      }
      
      for(int k = 0; k < 10; ++ k) local_top[k] = temp_top[k];
      
    }
    
    for(int k = 0; k < 10; ++ k) shared_top[tid][k] = local_top[k];
    __syncthreads();
    
    for(int d = nthreads / 2; d; d >>= 1){
      if(tid < d){
	int h1 = 0, h2 = 0, k;
	for(k = 0; k < 10; ++ k){
	  int d1 = shared_top[tid][h1], d2 = shared_top[tid + d][h2];
	  if(d1 == -1 && d2 == -1) break;
	  if(d1 > d2)
	    temp_top[k] = d1, ++ h1;
	  else
	    temp_top[k] = d2, ++ h2;
	}
	while(k < 10) temp_top[k] = -1, ++k;
	
	for(k = 0; k < 10; ++ k) shared_top[tid][k] = temp_top[k];
      }
      __syncthreads();
    }
    
    if(tid < 10)
      top[(v * n2 + o) * 10 + tid] = shared_top[0][tid];
  }
}

__global__ void init_data(int n1, int n2, int* top, int* weight_table){
  int tid = threadIdx.x, bid = blockIdx.x;

  for(int i = bid; i < n1; i += gridDim.x)
    for(int j = tid; j < n2; j += blockDim.x)
      top[(i * n2 + j) * 10] = weight_table[i * n2 + j];
}

//vector<int> I, O;
int main(int argc, char ** argv){
  
  int n, m, n1, n2;
  scanf("%d %d %d %d",&n,&m,&n1,&n2);

  int *node_weight = new int[n];
  int *node_type = new int[n];
  node **head = new node*[n];
  int *id = new int[n2 + 10];
  int *I = new int[n1 + 10];
  int cnt0 = 0, cnt1 = 0;
  for(int i = 0; i < n; ++ i){
    scanf("%d %d",node_weight + i, node_type + i);
    if(node_type[i] == 2) id[cnt1 ++] = i;
    if(node_type[i] == 1) I[cnt0 ++] = i;
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
  
  int *du = new int[n];
  int *edge_buf = new int[m];
  int *edge_offset = new int[n + 1];
  int *edge_inv_buf = new int[m];
  int *edge_inv_offset = new int[n + 1];
  int *sdu = new int[n];
  
  memset(du, 0, sizeof(du));

  double t0 = get_time();
  
  int offset = 0;
  for(int i = 0; i < n; ++ i){
    edge_offset[i] = offset;
    for(node *p = head[i]; p; p = p -> nxt){
      du[p -> p] ++;
      edge_buf[offset ++] = p -> p;
    }
  }

  edge_offset[n] = offset;
  for(int i = 0; i < n; ++ i) sdu[i] = (i ? sdu[i - 1] : 0) + du[i];
  for(int i = 0; i < n; ++ i) edge_inv_offset[i] = i ? sdu[i - 1] : 0;
  for(int i = 0; i < n; ++ i)
    for(node *p = head[i]; p; p = p -> nxt)
      edge_inv_buf[--sdu[p -> p]] = i;
  
  edge_inv_offset[n] = m;
  
  int* d_edge_offset;
  cudaMalloc(&d_edge_offset, sizeof(int) * (n + 1));
  int* d_edge_buf;
  cudaMalloc(&d_edge_buf, sizeof(int) * m);
  int* d_edge_inv_offset;
  cudaMalloc(&d_edge_inv_offset, sizeof(int) * (n + 1));
  int* d_edge_inv_buf;
  cudaMalloc(&d_edge_inv_buf, sizeof(int) * m);
  int* d_top;
  //cerr << sizeof(size_t) << endl;
  size_t sz = sizeof(int) * n * n2 * 10; 
  cerr << sz << endl;
  cudaMalloc(&d_top, sz);
  int* d_lock;
  cudaMalloc(&d_lock, sizeof(int) * n);
  int* d_inQue;
  cudaMalloc(&d_inQue, sizeof(int) * n);
  int* d_outQue;
  cudaMalloc(&d_outQue, sizeof(int) * n);
  int* d_du;
  cudaMalloc(&d_du, sizeof(int) * n);
  int* d_inQue_offset;
  cudaMalloc(&d_inQue_offset, sizeof(int) * n);
  int* d_sum;
  cudaMalloc(&d_sum, sizeof(int) * n);
  int* d_weight_table;
  cudaMalloc(&d_weight_table, sizeof(int) * n1 * n2);
  int* d_id;
  cudaMalloc(&d_id, sizeof(int) * n2);
  int* d_node_type;
  cudaMalloc(&d_node_type, sizeof(int) * n);
  int* d_node_weight;
  cudaMalloc(&d_node_weight, sizeof(int) * n);

  cudaMemcpy(d_edge_offset, edge_offset, sizeof(int) * (n + 1), cudaMemcpyHostToDevice);
  cudaMemcpy(d_edge_buf, edge_buf, sizeof(int) * m, cudaMemcpyHostToDevice);
  cudaMemcpy(d_edge_inv_offset, edge_inv_offset, sizeof(int) * (n + 1), cudaMemcpyHostToDevice);
  cudaMemcpy(d_edge_inv_buf, edge_inv_buf, sizeof(int) * m, cudaMemcpyHostToDevice);
  cudaMemcpy(d_du, du, sizeof(int) * n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_weight_table, weight_table, sizeof(int) * n1 * n2, cudaMemcpyHostToDevice);
  cudaMemcpy(d_id, id, sizeof(int) * n2, cudaMemcpyHostToDevice);
  cudaMemcpy(d_node_type, node_type, sizeof(int) * n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_node_weight, node_weight, sizeof(int) * n, cudaMemcpyHostToDevice);
  cudaMemset(d_top, -1, sz);
  cudaMemset(d_lock, -1, sizeof(int) * n);
  int in_size, out_size;

  in_size = n1;
  cudaMemset(d_inQue, 0, sizeof(int) * n);
  cudaMemset(d_outQue, 0, sizeof(int) * n);
  cudaMemcpy(d_inQue, I, sizeof(int) * n1, cudaMemcpyHostToDevice);
  
  //cudaMemcpy(d_inQue, I, sizeof(int) * n1, cudaMemcpyHostToDevice);
  //int inq[2010];
  //cudaMemcpy(inq, d_inQue, in_size * sizeof(int), cudaMemcpyDeviceToHost);
  //for(int i = 0 ; i < in_size; ++ i)
  //  cerr << inq[i] - I[i] << ' ';
  //cerr << endl;
  //int *lock = new int[n];

  /*
  cudaMemcpy(lock, d_lock, n * sizeof(int), cudaMemcpyDeviceToHost);
  for(int i = 0; i < n; ++ i)
    if(lock[i] != -1)
      cerr << i << ' ' << lock[i] << endl;
  */
  
  init_data<<<128,128>>>(n1, n2, d_top, d_weight_table); 
  int round = 0;
  
  while(true){
    cerr << "*** " << in_size << endl;
    expand<<<in_size, 32>>>(in_size, d_inQue, d_edge_offset, d_edge_buf, d_du, d_lock);
    count_inQue_offset<<<in_size, 32>>>(in_size, d_inQue, d_edge_offset, d_edge_buf, d_du, d_lock, d_inQue_offset);
    naive_scan<<<1,32>>>(in_size, d_inQue_offset, d_sum);
    
    cudaMemcpy(&out_size, d_sum + (in_size - 1), sizeof(int), cudaMemcpyDeviceToHost);
    enqueue<<<in_size,32>>>(in_size, d_inQue, d_outQue, d_edge_offset, d_edge_buf, d_lock, d_sum);
    cerr << "round " << ++ round << ": " << out_size << endl;
    if(out_size <= 0) break;
    update<<<256,nthreads>>>(out_size, n1, n2, d_node_weight, d_node_type, d_id, d_outQue, d_edge_inv_offset, d_edge_inv_buf, d_top);
    swap(in_size, out_size);
    swap(d_inQue, d_outQue);
  }

  double t1 = get_time();
  cerr << "runtime: " << t1 - t0 << endl;
  
  freopen("cuda_result.txt","w",stdout);
  int top[10];
  for(int i = 0; i < n2; ++ i){
    printf("%d ",id[i]);
    int x = id[i], cnt = 0;
    size_t offset = (1LL * x * n2 + i) * 10LL;
    cudaMemcpy(top, d_top + offset, sizeof(int) * 10, cudaMemcpyDeviceToHost);
    for(int j = 0; j < 10 && top[j] != -1; ++ j) ++ cnt;
    printf("%d ",cnt);
    for(int j = 0; j < cnt; ++ j)
      printf("%d ",top[j]);
    printf("\n");
  }
  
  return 0;
}

