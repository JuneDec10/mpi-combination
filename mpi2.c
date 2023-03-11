#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<mpi.h>
//#include <immintrin.h>

double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    double* rebuiltCoord = (double*)malloc(sizeof(double) * n * k);
    int i;
double chebyshevSum = 0;
for(int i=0; i<n; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            double distance = 0;
            int pivoti = pivots[ki];
            int j;
            for(j=0; j<dim; j++){
                distance += pow(coord[pivoti*dim + j] - coord[i*dim + j] ,2);
            
            }
                        
            rebuiltCoord[i*k + ki] = sqrt(distance);
        }
        
        for(int j=i; j>=0; j--){
            double chebyshev = 0;
            int ki;
            for(ki=0; ki<k; ki++){
                double dis = fabs(rebuiltCoord[i*k + ki] - rebuiltCoord[j*k + ki]);
                chebyshev = dis>chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
        
    }

    free(rebuiltCoord);

    return chebyshevSum;
}


void Combination(int ki, const int k, const int n, const int dim, const int M, double* coord, int* pivots,
                 double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){

     int start;
        int end;
        int size=-1;
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        int arrange = (n-pivots[ki-1]-1)/size;
        int yuarrange = (n-pivots[ki-1]-1)%size;

        int rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

       if(ki==k-1){
        int i;

      if(((pivots[ki-1])%2)==0){
       
        if(rank<yuarrange){
            start =(pivots[ki-1]+1)+ rank*(arrange+1);
            end = ((pivots[ki-1]+1)+((rank+1)*(arrange+1)) < n)?(pivots[ki-1]+1)+(((rank+1)*(arrange+1))):n;

        }else{

           start =(pivots[ki-1]+1)+  rank*arrange+yuarrange;
        end = ((pivots[ki-1]+1)+ ((rank+1)*arrange+yuarrange) < n)? (pivots[ki-1]+1)+ (((rank+1)*arrange+yuarrange)) : n ;
        }
        
      }else{
          
         if(rank>=(size-yuarrange)){
            start =(pivots[ki-1]+1)+ (size-rank-1)*(arrange+1);
            end = ((pivots[ki-1]+1)+(((size-rank-1)+1)*(arrange+1)) < n)?(pivots[ki-1]+1)+((((size-rank-1)+1)*(arrange+1))):n;
        }else{
           start =(pivots[ki-1]+1)+  (size-rank-1)*arrange+yuarrange;
        end = ((pivots[ki-1]+1)+ (((size-rank-1)+1)*arrange+yuarrange) < n)? (pivots[ki-1]+1)+ ((((size-rank-1)+1)*arrange+yuarrange)) : n ;
        }
      }

        for(i=start; i<end; i++){
            pivots[ki] = i;
            

            double distanceSum = SumDistance(k, n, dim, coord, pivots);
            
           
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            int kj;
            for(kj=0; kj<k; kj++){
                maxDisSumPivots[M*k + kj] = pivots[kj];
            }
            for(kj=0; kj<k; kj++){
                minDisSumPivots[M*k + kj] = pivots[kj];
            }

            int a;
            for(a=M; a>0; a--){
                if(maxDistanceSum[a] > maxDistanceSum[a-1]){
                    double temp = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a-1];
                    maxDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = maxDisSumPivots[a*k + kj];
                        maxDisSumPivots[a*k + kj] = maxDisSumPivots[(a-1)*k + kj];
                        maxDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
            //}
            //for(a=M; a>0; a--){
                if(minDistanceSum[a] < minDistanceSum[a-1]){
                    double temp = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a-1];
                    minDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = minDisSumPivots[a*k + kj];
                        minDisSumPivots[a*k + kj] = minDisSumPivots[(a-1)*k + kj];
                        minDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
            }
     
        }
        return;
    }

    int i;
    for(i=pivots[ki-1]+1; i<n; i++) {
        pivots[ki] = i;
        Combination(ki+1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);
    }
}

int main(int argc, char* argv[]){

    char* filename = (char*)"uniformvector-2dim-5h.txt";
    if( argc==2 ) {
        filename = argv[1];
    }  else if(argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }

    const int M = 1000;
    int dim;
    int n;
    int k;


    FILE* file = fopen(filename, "r");
    if( file == NULL ) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);


    struct timeval start;
    double* coord = (double*)malloc(sizeof(double) * dim * n);
    int i;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<dim; j++){
            fscanf(file, "%lf", &coord[i*dim + j]);
        }
    }
    fclose(file);
    gettimeofday(&start, NULL);


    double* maxDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    for(i=0; i<M; i++){
        maxDistanceSum[i] = 0;
    }

    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            maxDisSumPivots[i*k + ki] = 0;
        }
    }

    double* minDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    for(i=0; i<M; i++){
        minDistanceSum[i] = __DBL_MAX__;
    }

    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            minDisSumPivots[i*k + ki] = 0;
        }
    }

    int* temp = (int*)malloc(sizeof(int) * (k+1));
    temp[0] = -1;

    MPI_Init(&argc,&argv);

    int size=-1;
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    Combination(0, k, n, dim, M, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

        

    double * recMaxDistanceSum = (double*)malloc(sizeof(double) * (1001));
    double * resMaxDisSum = (double*)malloc(sizeof(double) * (1001));
    double * recMinDistanceSum = (double*)malloc(sizeof(double) * (1001));
    double * resMinDisSum = (double*)malloc(sizeof(double) * (1001));
    int * recMaxDistancePivots = (int*)malloc(sizeof(int) * (1001*k));
    int * resMaxPivots = (int*)malloc(sizeof(int) * (1001*k));
    int * recMinDistancePivots = (int*)malloc(sizeof(int) * (1001*k));
    int * resMinPivots = (int*)malloc(sizeof(int) * (1001*k));


int sizei = (size-1)/2;
int ri = size;
while(rank<ri){

    if(rank>sizei){ 
        MPI_Send(maxDistanceSum,1000,MPI_DOUBLE,rank-sizei-1,0,MPI_COMM_WORLD);
        MPI_Send(maxDisSumPivots,1000*k,MPI_INT,rank-sizei-1,1,MPI_COMM_WORLD);
        MPI_Send(minDistanceSum,1000,MPI_DOUBLE,rank-sizei-1,2,MPI_COMM_WORLD);
        MPI_Send(minDisSumPivots,1000*k,MPI_INT,rank-sizei-1,3,MPI_COMM_WORLD);
    }else{          
        if((rank==sizei)&&(ri%2==1)){            
        }else{
            MPI_Status status;
            MPI_Recv(recMaxDistanceSum,1000,MPI_DOUBLE,rank+sizei+1,0,MPI_COMM_WORLD,&status);
            MPI_Recv(recMaxDistancePivots,1000*k,MPI_INT,rank+sizei+1,1,MPI_COMM_WORLD,&status);
            MPI_Recv(recMinDistanceSum,1000,MPI_DOUBLE,rank+sizei+1,2,MPI_COMM_WORLD,&status);
            MPI_Recv(recMinDistancePivots,1000*k,MPI_INT,rank+sizei+1,3,MPI_COMM_WORLD,&status);
        
            int p1=0,p2=0;
            //bool ifp1=false,ifp2=false;
            for(int t = 0;t<1000;t++){
                if(maxDistanceSum[p1]>=recMaxDistanceSum[p2]){
                    resMaxDisSum[t] = maxDistanceSum[p1];
                    for(int g = 0;g<k;g++){
                    resMaxPivots[t*k+g]=maxDisSumPivots[p1*k+g];
                    }
                    p1++;
                    
                }else{
                    resMaxDisSum[t] = recMaxDistanceSum[p2];
                    for(int g = 0;g<k;g++){
                    resMaxPivots[t*k+g]=recMaxDistancePivots[p2*k+g];
                    }
                    p2++;
                }
            }
            
            for(int t = 0;t<1000;t++){
                maxDistanceSum[t] = resMaxDisSum[t];
                for(int g = 0;g<k;g++){
                    maxDisSumPivots[t*k+g]=resMaxPivots[t*k+g];
                }
            }

            p1=0;
            p2=0;
            for(int t = 0;t<1000;t++){
                if(minDistanceSum[p1]<=recMinDistanceSum[p2]){  
                    resMinDisSum[t] = minDistanceSum[p1];
                    for(int g = 0;g<k;g++){
                    resMinPivots[t*k+g]=minDisSumPivots[p1*k+g];
                    }
                    p1++;
                }else{
                    resMinDisSum[t] = recMinDistanceSum[p2];
                    for(int g = 0;g<k;g++){
                    resMinPivots[t*k+g]=recMinDistancePivots[p2*k+g];
                    }
                    p2++;
                }
            }
            
            for(int t = 0;t<1000;t++){
                minDistanceSum[t] = resMinDisSum[t];
                for(int g = 0;g<k;g++){
                	maxDisSumPivots[t*k+g]=resMaxPivots[t*k+g];
                    minDisSumPivots[t*k+g]=resMinPivots[t*k+g];
                }
            }
        }
    
    }
    ri = sizei+1;
    if(sizei==0) break;
    sizei = sizei/2; 
};
 
    struct timeval end;
    gettimeofday (&end, NULL);
    printf("rank %d Using time : %f ms\n", rank,(end.tv_sec-start.tv_sec)*1000.0+(end.tv_usec-start.tv_usec)/1000.0);

    if(rank==0){
        FILE* out = fopen("result.txt", "w");
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", maxDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i*k + k-1]);
    }
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", minDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i*k + k-1]);
    }
    fclose(out);
    }
    

    MPI_Finalize();
    return 0;
}
