%%cuda --name r3_parallel.cu

// We will only be using static arrays in this
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

const double INF_ = 1e9;
const double epsilon = 0.0000000001;

struct ii {
    int F;
    int S;
};

__device__ double calc_dist(ii p1, ii p2){
    double x1 = p1.F; double y1 = p1.S;
    double x2 = p2.F; double y2 = p2.S;
    double dx = x2-x1; double dy = y2-y1;
    return sqrt((dx*dx) + (dy*dy));
}


__device__ bool isvalid(ii p, int n, int m, double dist_from_origin, double range){
    if(p.F < 0 || p.S < 0) return 0;
    if(p.F >= n || p.S >= m) return 0;
    if(dist_from_origin > range) return 0;
    return 1;
}


__device__ double find_slope(int *grid, ii observer, ii curr, int observer_ht, int n, int m){
    double dist = calc_dist(observer, curr);
    int id_curr = (curr.F)*(m) + (curr.S);
    int id_obs = (observer.F)*(m) + (observer.S);
    double dh = double(grid[id_curr]) - (double(grid[id_obs]) + (double)observer_ht);
    double m_curr = dh/dist;
    return m_curr;
}


__device__  void plotPixel(ii *line, int &numpts, int x1, int y1, int x2, int y2, int dx, int dy, int decide, int lid_add){
    int pk = 2 * dy - dx;
    numpts = 0; 

    for (int i = 0; i <= dx; i++) {
        x1 < x2 ? x1++ : x1--;
        if (pk < 0) {
            if (decide == 0) {
                ii temp = {x1,y1};
                line[lid_add+numpts] = temp;
                numpts++;
                pk = pk + 2 * dy;
            }
            else {
                ii temp = {y1,x1};
                line[lid_add+numpts] = temp;
                numpts++;
                pk = pk + 2 * dy;
            }
        }
        else {
            y1 < y2 ? y1++ : y1--;
            if (decide == 0) {
                ii temp = {x1,y1};
                line[lid_add+numpts] = temp;
                numpts++;
            }
            else {
                ii temp = {y1,x1};
                line[lid_add+numpts] = temp;
                numpts++;
            }
            pk = pk + 2 * dy - 2 * dx;
        }
    }
}


__device__ void bresenham(ii *line_grid, int &numpts, int x1, int y1, int x2, int y2, int lid_add){
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    numpts=0;
    if(dx > dy){
        plotPixel(line_grid, numpts, x1, y1, x2, y2, dx, dy, 0, lid_add);
    }
    else{
        plotPixel(line_grid, numpts, y1, x1, y2, x2, dy, dx, 1, lid_add);
    }
}


void readfile_and_make_grid(int* grid, int &n, int &m){
    //std::ifstream infile("/Users/avirajbevli/Desktop/Sem10/HPPP/cuda_term project/data/cdnf43w/heights.txt");
    std::cout<<"Reading from input file..."<<std::endl;
    std::ifstream infile("heights.txt");
    if(!infile.is_open()) {
        printf("not abel to read \n");
        std::cerr << "Error: Could not open read file" << std::endl;
        return  ;
    }

    std::string line;
    n=0;
    int numpts=0;
    while(std::getline(infile, line)) {
        std::string str; 
        int i=0;
        while(i<line.size()){
            if(line[i]==','){
                int temp = stoi(str);
                grid[numpts++] = temp;
                str="";
                i++;
            }
            else{
                str+=line[i];
            }
            i++;
        }
        n++;
    }
    m = (numpts)/n;
    infile.close();
    std::cout<<"Done reading \n"<<std::endl;
}


// cant rely on print in global
__global__ void test_kernel(int *viewshed_d, int *grid_d){
    viewshed_d[0]=999;
    grid_d[0]=999;
}

// -1: out of range, 0: not visible, 2: visible
// can not print anything inside global/device functions
__global__ void calculate_viewshed_and_writefile_r3(int *grid, int* viewshed, ii *line_grid, int n, int m, int range){    
    // viewshed[0]=999;
    // grid[0]=999;

    int obsy = 1000; int obsx = 3000;
    ii observer = {obsy, obsx}; //<y,x> or <row,col>
    int observer_ht = 10; // height of observer from the ground at <observer_x, observer_y>  

    int I = observer.F + (blockIdx.y * blockDim.y + threadIdx.y - range); 
    int J = observer.S + (blockIdx.x * blockDim.x + threadIdx.x - range);
    if(I>=0 && J>=0 && I<n && J<m){
        ii target = {I,J};
        int index = (target.F)*(m) + (target.S);
        viewshed[index]=2;
        __syncthreads(); // no need for this, right?
        
        double m_target = find_slope(grid, observer, target, observer_ht, n, m); 

        int linesz=0;
        int tempy = blockIdx.y * blockDim.y + threadIdx.y;
        int tempx = blockIdx.x * blockDim.x + threadIdx.x;
        int lid_add = (tempy*(2*range + 1) + tempx)*(2*range); // range
        bresenham(line_grid, linesz, observer.F, observer.S, target.F, target.S, lid_add);

        // printf("I:%d, J:%d, linesz:%d \n", I, J, linesz);
        double m_max = -INF_;
        for(int lid=1;lid<linesz-1;lid++){
            double m_curr = find_slope(grid, observer, line_grid[lid_add+lid], observer_ht, n, m);
            m_max = max(m_max, m_curr);
        }
        if(m_max > m_target + epsilon){
            viewshed[index]=0;
        }
        __syncthreads(); // so that the viewshed is complete
    }
}


__global__ void kernel_r3(int *grid, int* viewshed, ii *line_grid, int n, int m, int range){    
    int obsy = 1000; int obsx = 3000;
    ii observer = {obsy, obsx}; //<y,x> or <row,col>
    int observer_ht = 10; // height of observer from the ground at <observer_x, observer_y>  

    int I = observer.F + (blockIdx.y * blockDim.y + threadIdx.y - range); 
    int J = observer.S + (blockIdx.x * blockDim.x + threadIdx.x - range);
    if(I>=0 && J>=0 && I<n && J<m){
        ii target = {I,J};
        int index = (target.F)*(m) + (target.S);
        viewshed[index]=1;
        __syncthreads(); // no need for this, right?
        
        double m_target = find_slope(grid, observer, target, observer_ht, n, m); 

        int linesz=0;
        int tempy = blockIdx.y * blockDim.y + threadIdx.y;
        int tempx = blockIdx.x * blockDim.x + threadIdx.x;
        int lid_add = (tempy*(2*range + 1) + tempx)*(2*range); // range
        bresenham(line_grid, linesz, observer.F, observer.S, target.F, target.S, lid_add);

        // printf("I:%d, J:%d, linesz:%d \n", I, J, linesz);
        double m_max = -INF_;
        for(int lid=1;lid<linesz-1;lid++){
            double m_curr = find_slope(grid, observer, line_grid[lid_add+lid], observer_ht, n, m);
            m_max = max(m_max, m_curr);
        }
        if(m_max > m_target + epsilon){
            if(viewshed[index]!=2){
                viewshed[index]=0;
            }
        }
        else{
            viewshed[index]=2;
        }
        __syncthreads(); // so that the viewshed is complete
    }
}


__global__ void kernel_r3_thread_coarsening(int *grid, int* viewshed, ii *line_grid, int n, int m, int range){    
    int obsy = 1000; int obsx = 3000;
    ii observer = {obsy, obsx}; //<y,x> or <row,col>
    int observer_ht = 10; // height of observer from the ground at <observer_x, observer_y>  

    int Ivec[4];
    int Jvec[4];
    Ivec[0] = observer.F + (blockIdx.y * 2 * blockDim.y + threadIdx.y - range); 
    Jvec[0] = observer.S + (blockIdx.x * 2 * blockDim.x + threadIdx.x - range);
    Ivec[1] = (Ivec[0] + blockDim.y); 
    Jvec[1] = Jvec[0];
    Ivec[2] = Ivec[0]; 
    Jvec[2] = (Jvec[0] + blockDim.x);
    Ivec[3] = (Ivec[0] + blockDim.y); 
    Jvec[3] = (Jvec[0] + blockDim.x);

    for(int ctrid=0;ctrid<4;ctrid++){
    	int I = Ivec[ctrid];
    	int J = Jvec[ctrid];

    	if(I>=0 && J>=0 && I<n && J<m){
	        ii target = {I,J};
	        int index = (target.F)*(m) + (target.S);
	        viewshed[index]=1;
	        __syncthreads(); // no need for this, right?
	        
	        double m_target = find_slope(grid, observer, target, observer_ht, n, m); 

	        int linesz=0;
	        int tempy = blockIdx.y * blockDim.y + threadIdx.y;
	        int tempx = blockIdx.x * blockDim.x + threadIdx.x;
	        int lid_add = (tempy*(2*range + 1) + tempx)*(2*range); // range
	        bresenham(line_grid, linesz, observer.F, observer.S, target.F, target.S, lid_add);

	        // printf("I:%d, J:%d, linesz:%d \n", I, J, linesz);
	        double m_max = -INF_;
	        for(int lid=1;lid<linesz-1;lid++){
	            double m_curr = find_slope(grid, observer, line_grid[lid_add+lid], observer_ht, n, m);
	            m_max = max(m_max, m_curr);
	        }
	        if(m_max > m_target + epsilon){
	            if(viewshed[index]!=2){
	                viewshed[index]=0;
	            }
	        }
	        else{
	            viewshed[index]=2;
	        }
	        __syncthreads(); // so that the viewshed is complete
    	}
    }
}


__global__ void kernel_r3_block_coarsening(int *grid, int* viewshed, ii *line_grid, int n, int m, int range){    
	int obsy = 1000; int obsx = 3000;
    ii observer = {obsy, obsx}; //<y,x> or <row,col>
    int observer_ht = 10; // height of observer from the ground at <observer_x, observer_y>  

    int Ivec[4];
    int Jvec[4];
    Ivec[0] = observer.F + (blockIdx.y * blockDim.y + threadIdx.y - range); 
    Jvec[0] = observer.S + (blockIdx.x * blockDim.x + threadIdx.x - range);
    Ivec[1] = (Ivec[0] + (gridDim.y)*(blockDim.y)); 
    Jvec[1] = Jvec[0];
    Ivec[2] = Ivec[0]; 
    Jvec[2] = (Jvec[0] + (gridDim.x)*(blockDim.x));
    Ivec[3] = (Ivec[0] + (gridDim.y)*(blockDim.y)); 
    Jvec[3] = (Jvec[0] + (gridDim.x)*(blockDim.x));

    for(int ctrid=0;ctrid<4;ctrid++){
    	int I = Ivec[ctrid];
    	int J = Jvec[ctrid];

    	if(I>=0 && J>=0 && I<n && J<m){
	        ii target = {I,J};
	        int index = (target.F)*(m) + (target.S);
	        viewshed[index]=1;
	        __syncthreads(); // no need for this, right?
	        
	        double m_target = find_slope(grid, observer, target, observer_ht, n, m); 

	        int linesz=0;
	        int tempy = blockIdx.y * blockDim.y + threadIdx.y;
	        int tempx = blockIdx.x * blockDim.x + threadIdx.x;
	        int lid_add = (tempy*(2*range + 1) + tempx)*(2*range); // range
	        bresenham(line_grid, linesz, observer.F, observer.S, target.F, target.S, lid_add);

	        // printf("I:%d, J:%d, linesz:%d \n", I, J, linesz);
	        double m_max = -INF_;
	        for(int lid=1;lid<linesz-1;lid++){
	            double m_curr = find_slope(grid, observer, line_grid[lid_add+lid], observer_ht, n, m);
	            m_max = max(m_max, m_curr);
	        }
	        if(m_max > m_target + epsilon){
	            if(viewshed[index]!=2){
	                viewshed[index]=0;
	            }
	        }
	        else{
	            viewshed[index]=2;
	        }
	        __syncthreads(); // so that the viewshed is complete
    	}
    }
}

int main(){
    cudaError_t err = cudaSuccess;
    int n=3600; int m=3600;
    int range = 200;
    size_t size = (3600*3600)*sizeof(int);
    // allocating a 2 dim kernel sizes 
    
    int gridCols = ceil(((2*200)+1)/16);
    int gridRows = ceil(((2*200)+1)/16);
    std::cout<<"gridCols: "<<gridCols<<", gridRows: "<<gridRows<<std::endl;
    dim3 gridDim(gridCols, gridRows, 1);
    dim3 blockDim(16, 16, 1); 

    int *grid_h = (int *)malloc(size);
    if (grid_h == NULL) { // check if memory allocation was successful
        printf("Memory allocation to grid failed.");
        return 1;
    }
    readfile_and_make_grid(grid_h,n,m);
    printf("read input file grid_0: %d \n", grid_h[0]);

    int *grid_d = NULL;
    cudaMalloc((void **) &grid_d, size);
    if (err != cudaSuccess){
        printf("FAIL");
        fprintf(stderr, "Failed to allocate device grid (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
  
    int *viewshed_h = (int *)malloc(size);
    for(int i =0;i<n;i++){
        for(int j=0;j<m;j++){
            viewshed_h[i*m+j] = -1;
        }
    }

    int *viewshed_d = NULL;
    cudaMalloc((void **) &viewshed_d, size);
    if (err != cudaSuccess){
        printf("FAIL");
        fprintf(stderr, "Failed to allocate device viewshed (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    ii *line_grid = NULL;
    size_t size_line_grid = (2*range + 1)*(2*range + 1)*(2*range)*sizeof(ii);
    cudaMalloc((void **) &line_grid, size_line_grid);
    if (err != cudaSuccess){
        printf("FAIL");
        fprintf(stderr, "Failed to allocate device viewshed (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(grid_d, grid_h, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to copy grid from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(viewshed_d, viewshed_h, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to copy viewshed from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    printf("before kernel call, viewshed_0: %d, grid_0: %d \n", viewshed_h[0], grid_h[0]);
    kernel_r3<<<gridDim, blockDim>>>(grid_d, viewshed_d, line_grid, n, m, range);
    err = cudaGetLastError();
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to launch kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    dim3 gridDim_thread_coarsening(gridCols, gridRows, 1);
    dim3 blockDim_thread_coarsening(8, 8, 1); 
    kernel_r3_thread_coarsening<<<gridDim_thread_coarsening, blockDim_thread_coarsening>>>(grid_d, viewshed_d, line_grid, n, m, range);
    err = cudaGetLastError();
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to launch kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    dim3 gridDim_block_coarsening(gridCols/2, gridRows/2, 1);
    dim3 blockDim_block_coarsening(16, 16, 1); 
    kernel_r3_block_coarsening<<<gridDim_block_coarsening, blockDim_block_coarsening>>>(grid_d, viewshed_d, line_grid, n, m, range);
    err = cudaGetLastError();
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to launch kernel (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(viewshed_h, viewshed_d, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to copy viewshed from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    printf("after kernel call, viewshed_0: %d, grid_0: %d \n", viewshed_h[0], grid_h[0]);


    std::cout<<"Writing to output viewshed file ...."<<std::endl;
    std::ofstream outFile("viewshed_r3.txt");   
    if (!outFile.is_open()) {
        std::cout << "Unable to open write file" << std::endl;
        return;
    }

    for(int i=0;i<n*m;i++){
        outFile<<viewshed_h[i]<<" ";
        if((i%m)==m-1){
            outFile<<std::endl;
        }
    }
    // Close the file
    outFile.close();
    std::cout<<"Done writing"<<std::endl;  

    // free device memory
    err = cudaFree(grid_d); 
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to free device grid (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaFree(viewshed_d); 
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to free device grid (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // reset the device and exit
    err = cudaDeviceReset();
    if (err != cudaSuccess){
        fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    //free host memory
    free(grid_h);
    free(viewshed_h);
    printf("Done \n");
    return 0;
}
