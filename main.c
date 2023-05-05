#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;

extern double func  (double t, double x);
extern double fi    (double x);
extern double ksi   (double t);
extern void   shema(int k, int m, double** u);

void print_time(int size, double t1) {
	double t2 = MPI_Wtime();
    	printf("%d %lg\n", size, t2 - t1);
    	FILE *fp;
   	fp = fopen("time.csv", "a");
    	fprintf(fp, "%d,%lg\n", size, t2 - t1);
    	fclose(fp);
}

int main(int argc, char** argv) {

    	int M = x_max / x_step; 
    	int K = t_max / t_step; 

    	double t1, t2;   
	int rank, size;    

	double **u = (double **)calloc(K, sizeof(double *));
        for (int i = 0; i < K; i++)
                u[i] = (double *)calloc(M, sizeof(double));

	MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    	MPI_Comm_size( MPI_COMM_WORLD, &size );
    
    	if(rank == 0)
    		t1 = MPI_Wtime();      

    	for (int m = 0; m < M; ++m)
        	u[0][m] = fi(m);         

    	for (int k = 0; k < K; ++k)
        	u[k][0] = ksi(k);        

	if (size == 1)
	{      
                for (int k = 0; k < K - 1; k++) {
                        for (int m = 1; m < M; m++)
                                shema(k, m, u);
                }

                print_time(size, t1);
                printf("result is in the file time.csv\n");
                        FILE *fp = fopen("output.csv", "w");
                        fprintf(fp, "x,t,u\n");
                        for (int k = 0; k < K; k++) {
                                double t = k * t_step;
                                for (int m = 0; m < M; m++) {
                                        double x = m * x_step;
                                        fprintf(fp, "%.6f,%.6f,%.6f\n", x, t, u[k][m]);
                                }
                        }
                fclose(fp);
        }

	else{  
        	int part = M / size; 
        	int shift = M % size;   
		
		int num_knots;
		int x_0;
		if (rank == 0){
			num_knots = part + shift;
			x_0 = 0;
		}
		else{
			num_knots = part;
			x_0 = part * rank + shift;
		}
        	
		int x_1 = x_0 + num_knots;       


        	for (int k = 0; k < K - 1; k++){
        		for (int m = x_0 + 1; m < x_1 - 1; m++) { 
        			shema(k, m, u);            
        		}
        		if (rank != 0 && rank != size - 1) {
        			if (k > 0){
        				MPI_Recv(&u[k][x_0 - 1], 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
    		    			MPI_Recv(&u[k][x_1], 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        			}
    	    			shema(k, x_0, u);     
    	    			shema(k, x_1 - 1, u);
    	    			if (k < K - 2){
    	    				MPI_Send(&u[k + 1][x_0], 1, MPI_DOUBLE, rank - 1, k + 1, MPI_COMM_WORLD);      
    		    			MPI_Send(&u[k + 1][x_1 - 1], 1, MPI_DOUBLE, rank + 1, k + 1, MPI_COMM_WORLD);
    	    			}
        		}	
        		else if (rank == 0) {
        			if (k > 0)
    	    				MPI_Recv(&u[k][x_1], 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      

    	    			shema(k, x_1 - 1, u);

    	    			if (k < K - 2)
    		    			MPI_Send(&u[k + 1][x_1 - 1], 1, MPI_DOUBLE, rank + 1, k + 1, MPI_COMM_WORLD);
        		}
        		else { 
        			if (k > 0)
    	    				MPI_Recv(&u[k][x_0 - 1], 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   

    	    			shema(k, x_0, u);

    	    			if(k < K - 2)
    		    			MPI_Send(&u[k + 1][x_0], 1, MPI_DOUBLE, rank - 1, k + 1, MPI_COMM_WORLD);
        		}
        		MPI_Barrier(MPI_COMM_WORLD); 
        	}

        	if (rank == 0) {  // окончательные данные
        		for (int r = 1; r < size; r++) {
                		for (int k = 1; k < K; k++){
                    			int first = part * r + shift;
                    			int last = first + part;
                    			MPI_Recv(&u[k][first], part, MPI_DOUBLE, r, K, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                		}
            		}
        	}

        	else { 
        		for (int k = 1; k < K; k++) 
        		 	MPI_Send(&u[k][x_0], num_knots, MPI_DOUBLE, 0, K, MPI_COMM_WORLD);    
        	}

        	if (rank == 0) {
            		print_time(size, t1);
	    	
                	printf("result is in the file time.csv\n");
                	FILE *fp = fopen("output.csv", "w");
                	fprintf(fp, "x,t,u\n"); 
                	for (int k = 0; k < K; k++) {
                    		double t = k * t_step;
                    		for (int m = 0; m < M; m++) {
                        		double x = m * x_step;
                        		fprintf(fp, "%.6f,%.6f,%.6f\n", x, t, u[k][m]);
                    		}
                	}
                fclose(fp);
            	}

    	} 
    		
	for (int i = 0; i < K; i++) 
        	free(u[i]);
    	free(u);

    	MPI_Finalize();
    
    	return 0;
}
