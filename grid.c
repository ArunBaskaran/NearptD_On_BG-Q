#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include <math.h>
//#include<pthreads.h>
// nx, ny, and nz are the nodes of the grid 

#define NX 10
#define NY 10
#define NZ 10
#define dx -1.0
#define dy -1.0
#define dz -1.0
#define FAILURE_VALUE 0.0
float queryX, queryY;

// A data structure to define cells that together form the grid.

struct Cell
{
	int xbegin, xend, ybegin, yend, zbegin, zend ;
	int * p_ids ;  // A dynamic array to store the ids of points in the cell
	int p_flag ; // An iterative variable that will give the number of points in the cell 
} ;


// A data structure to store 3d points. The points will be stored in an array of objets of this structure, say Point_3D point[points_max], and the array index will be the point's integer id.


struct Point_3D
{
	double x, y, z ;
} ;

//Function declarations
double distance(struct Point_3D point1, struct Point_3D point2) ;

struct Neighbor
{
	float distance;
	int id;
};

int i, j, k ;
int Lx, Ly, Lz ;
int cell_flag ;
int mpi_myrank;
int mpi_commsize;
int cell_max = NX * NY * NZ ;
int numPoints ;
int numPoints_rank ;
int myChunkStart ;
int myChunkEnd ;
int numQueryPoints ;
int points_in_cell[NX * NY * NZ] ;
struct Point_3D * allPoints ;
struct Cell * cells ;


MPI_File fh_input ;
MPI_File fh_query ;
MPI_File fh_output ;

char* inFileName ;
char* queryFileName ;
char* outFileName ;

void grid_init()
{
	cell_flag = mpi_myrank * (cell_max/mpi_commsize) ; 
	for(i = mpi_myrank*(NX/mpi_commsize) ; i < (mpi_myrank+1)*(NX/mpi_commsize) ; i++)
	{
		for(j = mpi_myrank*(NY/mpi_commsize) ; j < (mpi_myrank+1)*(NY/mpi_commsize) ; j++)
		{
			for(k = mpi_myrank*(NZ/mpi_commsize) ; k < (mpi_myrank+1)*(NZ/mpi_commsize); k++)
			{
				cells[cell_flag].xbegin = i*dx ;
				cells[cell_flag].xend = (i+1)*dx ;
				cells[cell_flag].ybegin = j*dy ;
				cells[cell_flag].yend = (j+1)*dy ;
				cells[cell_flag].zbegin = k*dz ;
				cells[cell_flag].zend = (k+1)*dz ;
				cells[cell_flag].p_ids = malloc(1*sizeof(struct Cell)) ;
				cell_flag++ ;
			}
		}
	}
}

void readInputFile()
{
	MPI_Offset offset = (int)sizeof(double)*mpi_myrank*numPoints_rank*3 ;
	double x,y,z; 
	for(i = 0; i < numPoints_rank; i++)
  	{
		offset += (int)sizeof(double)*i*3  ;
		MPI_File_read_at_all(fh_input , offset , &x , 1 , MPI_FLOAT , MPI_STATUS_IGNORE);
		MPI_File_read_at_all(fh_input , offset , &y , 1 , MPI_FLOAT , MPI_STATUS_IGNORE);
		MPI_File_read_at_all(fh_input , offset , &z , 1 , MPI_FLOAT , MPI_STATUS_IGNORE);		
  		//fscanf(inFile, "%f", &x);
  		//fscanf(inFile, "%f", &y);
  		//fscanf(inFile, "%f", &z);

		//Mutex lock ??
  		allPoints[mpi_myrank*numPoints_rank + i].x = x;
  		allPoints[mpi_myrank*numPoints_rank + i].y = y;
  		allPoints[mpi_myrank*numPoints_rank + i].z = z;        
  	}
}

//A function that identifies a point's cell, and updates the point's id to the cell's list of points. 
void home_cell(int id)
{ 
	int i ;
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;
	for(int i = 0 ; i < cell_max ; i++)
	{
		if((x1>(double)cells[i].xbegin) & (x1<(double)cells[i].xend) & (y1>(double)cells[i].ybegin) & (y1<(double)cells[i].yend) & (z1>(double)cells[i].zbegin) & (z1<(double)cells[i].zend))
		{
			if(cells[i].p_flag==0)
			{
				//cells[i].p_ids = malloc(1*sizeof(struct Cell)) ;
				cells[i].p_ids[0] = id ;
				points_in_cell[i] = 1 ;
				//exit(0) ;
			}
			else
			{		
				cells[i].p_ids = realloc(cells[i].p_ids, 1*sizeof(struct Cell)) ;
				cells[i].p_ids[points_in_cell[i]] = id ;
				points_in_cell[i]++ ;
				//exit(0) ;
			}
		}
	}
}

int get_cell(int id)
{
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;
	for(int i = 0 ; i < cell_max ; i++)
	{
		if((x1>(double)cells[i].xbegin) & (x1<(double)cells[i].xend) & (y1>(double)cells[i].ybegin) & (y1<(double)cells[i].yend) & (z1>(double)cells[i].zbegin) & (z1<(double)cells[i].zend))
		{
			return i ;
		}
	}
}


void assignPointsToCells()
{
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;

  	for(int i = mpi_myrank * numPoints_rank ; i < (mpi_myrank+1) * numPoints_rank ; i++)
  	{
  		home_cell(i);
	}
}	

//Function to read from query input file. Stores in queryX, queryY
//Only rank = 0 reads this
//Input file should be space separated in input/input.txt
void readQueryFile()
{
	if (mpi_myrank == 0)
	{
		FILE *fp;
		fp = fopen("input/input.txt","r");
		fscanf(fp,"%f",&queryX);
		fscanf(fp,"%f",&queryY);
	}
}


double distance(struct Point_3D point1, struct Point_3D point2)
{

	double distance=0.0;

	distance= sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2));


	return distance;
}

//TODO add nearest neighbor logic (3d modular arithmetic.)

//Assuming that index is the index of the query point

struct Neighbor NearestNeighbor(int index)
{
	//struct Neighbor n = malloc(sizeof(Neighbor));
	double min_distance = 100000 ;
	int home_cell_id = get_cell(index) ;
	int index_new, neighbor_id ;
	int success_flag = 0 ;
	double temp_distance = 100000 ;
	for(int i = -1 ; i<=1 ; i++)
	{
		for(int j = -1 ; j <= 1 ; j++)
		{
			for(int k = -1 ; k <= 1 ; k++)
			{
				if(i==0 & j==0 & k==0)
					{continue ; }
				else
				{ 
					index_new = index + i + j*NY + k*(NX*NY) ; 
					if((index_new >= 0) & (index_new < cell_max))
					{ 
						int temp = 0;
						while(temp <= points_in_cell[index_new])
						{					

							if(points_in_cell[index_new]>0)
							{    
								int temp_id = cells[index_new].p_ids[temp] ;
							
								temp_distance = distance(allPoints[index], allPoints[temp_id]) ;
								if(temp_distance < min_distance)
								{
									min_distance = temp_distance ;
									neighbor_id = cells[i].p_ids[temp] ;
									success_flag++ ;
								}
							}
							temp++ ; 
						}
					}
				}

			}
		}
	}



	// Neighboring cells are : index+1, index+1-(NX*NY), index+1+(NX*NY), index+1+NY, index+1-NY, index+1-NY-(NX*NY), index+1-NY+(NX*NY), index+1+NY-(NX*NY), index+1+NY-(NX*NY),///// index-1, index-1-(NX*NY), index-1+(NX*NY), index-1+NY, index-1-NY, index-1-NY-(NX*NY), index-1-NY+(NX*NY), index-1+NY-(NX*NY), index-1+NY-(NX*NY), index+NY, index-NY, index + (NX * NY) , index - (NX * NY) //// index-NY-(NX*NY), index-NY+(NX*NY), index+NY-(NX*NY), index+NY-(NX*NY)   (26 cells for 3D)

	struct Neighbor n ;
	if(success_flag==0)
	{
		n.id = 0.0 ;
		n.distance = FAILURE_VALUE ;
	}
	else
	{
		n.id = neighbor_id ;
		n.distance = min_distance ;
	}
	return n;	
}


struct Neighbor NearestNeighborExhaustive(int index, int myChunkStart, int myChunkEnd)
{
	float distance_final = 10000000;
	int nearestNeighbor = FAILURE_VALUE;
	for(int i = myChunkStart; i < myChunkEnd; i++)
	{
		if(index == i)
		{
			continue;
		}
		float temp = distance(allPoints[index], allPoints[i]);
		if(temp < distance_final)
		{
			distance_final = temp;
			nearestNeighbor = i;
		}
	}
	struct Neighbor n ;
	n.distance = distance_final;
	n.id = nearestNeighbor;
	return n;
}

//args
//1: total # input points
//2: inputFile.txt
//3: total # input query points
//4: queryFile.txt
//5: outFile.txt
int main(int argc, char** argv)
{
	MPI_Comm file_comm ;
	Lx = dx * NX ;
	Ly = dy * NY ;
	Lz = dz * NZ ;
	cell_flag = 0 ;

	if(argc < 6)
	{
		perror("Missing arguments. Exepects 1) inputFile.txt ") ;
	}
	numPoints = atoi(argv[1]);
	inFileName = argv[2];
	numQueryPoints = atoi(argv[3]);
	queryFileName = argv[4];
	outFileName = argv[5] ;

        //MPI_File_open(file_comm,inFileName, MPI_MODE_CREATE | MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_input);

	cells = malloc((NX*NY*NZ) * sizeof(struct Cell)) ;
	allPoints = malloc(numPoints * sizeof(struct Point_3D)) ;

  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;
    
        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

	numPoints_rank = numPoints/mpi_commsize ;
	grid_init() ;
	for(int i = 0 ; i < cell_max ; i++)
	{
		cells[i].p_flag = 0 ;
	}


        double startTime;
        double endTime;
        if(mpi_myrank == 0)
        {
    		startTime = MPI_Wtime(); //rank 0 is the timekeeper.
        }

	readInputFile() ;
	assignPointsToCells() ;

//	FILE *inFile = fopen(inFileName, "r");
	FILE *queryFile = fopen(queryFileName, "r");
	FILE *outFile = fopen(outFileName, "w");


        //point = malloc(sizeof(point3D) * numPoints);

  	//int chunkSize = (numPoints / mpi_commsize);



	int* queryPoints = malloc(sizeof(int) * numQueryPoints);
	for(i = 0; i < numQueryPoints; i++)
  	{
  		int index;
  		fscanf(queryFile, "%d", &index);
  		queryPoints[i] = index;
  	}

	int numQueryPoints_rank = numQueryPoints/mpi_commsize ;
//To-do : Refine the query points
  	for(int i = mpi_myrank*numQueryPoints_rank ; i < (mpi_myrank+1)*numQueryPoints_rank ; i++) //for every query point (basic version 2 search modes.)
	//Should this loop over numQueryPoints/mpi_commsize number of points?
  	{
  		struct Neighbor n = NearestNeighbor(queryPoints[i]);  
  		float minDistance;
  		minDistance = n.distance;
  		//MPI_Allreduce(&minDistance, &minDistance, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  		/*if(n.distance == (float)FAILURE_VALUE)
  		{
  			struct Neighbor n = NearestNeighborExhaustive(queryPoints[i], myChunkStart, myChunkEnd);
  			minDistance = n.distance;
			MPI_Barrier(MPI_COMM_WORLD) ;
  			MPI_Allreduce(&minDistance, &minDistance, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

  			if(minDistance == n.distance)
  			{
  				//TODO we need some kind of locking for ties.
  				if(n.distance == (float)FAILURE_VALUE)
  				{
  					if(mpi_myrank == 0)
  					{
  						fprintf(outFile , "%d nearest neighbor: Nonexistent\n",queryPoints[i]);
  					}
  				}
  				else
  				{
  					fprintf(outFile , "%d nearest neighbor: %d\n",queryPoints[i], n.id) ;
  				}
  			}
  		}*/
  	}

	if(mpi_myrank == 0)
    {
    	endTime = MPI_Wtime();
    	double totalTime = endTime-startTime;
    	printf("Total elapsed time was %lf\n", totalTime);
    } 


	MPI_File_close(&fh_input) ;
	MPI_File_close(&fh_query) ;
	MPI_File_close(&fh_output) ;
	free(cells) ;
	free(allPoints) ;

  	MPI_Finalize();

}




  //Added by Sidharth

		
