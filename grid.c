#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

// nx, ny, and nz are the nodes of the grid 

#define NX 10
#define NY 10
#define NZ 10
#define dx 0.1
#define dy 0.1
#define dz 0.1
#define FAILURE_VALUE = -1;

// A data structure to define cells that together form the grid.

struct Cell
{
	int xbegin, xend, ybegin, yend, zbegin, zend ;
	int * p_ids ;  // A dynamic array to store the ids of points in the cell
	int p_flag = 0 ; // An iterative variable that will give the number of points in the cell 
}


// A data structure to store 3d points. The points will be stored in an array of objets of this structure, say Point_3D point[points_max], and the array index will be the point's integer id.

struct Point_3D
{
	double x, y, z ;
}

struct Neighbor
{
	float distance;
	int id;
};

int i, j, k ;
int Lx, Ly, Lz ;
Lx = dx * NX ;
Ly = dy * NY ;
Lz = dz * NZ ;
int cell_flag = 0 ;
const int cell_max = NX * NY * NZ ;

Point_3D point*; //A list of all points
Cell cells[cell_max] ;

for(i = 0 ; i < NX ; i++)
{
	for(j = 0 ; j < NY ; j++)
	{
		for(k = 0 ; k < NZ ; k++)
		{
			cells[cell_flag].xbegin = i ;
			cells[cell_flag].xend = i+1 ;
			cells[cell_flag].ybegin = j ;
			cells[cell_flag].yend = j+1 ;
			cells[cell_flag].zbegin = k ;
			cells[cell_flag].zend = k+1 ;
			cell_flag++ ;
		}
	}
}

//A function that identifies a point's cell, and updates the point's id to the cell's list of points. 
void home_cell(int id)
{ 
	int i ;
	double x1 = point[id].x ;
	double y1 = point[id].y ;
	double z1 = point[id].z ;
	for(int i = 0 ; i < cell_max ; i++)
	{
		if((x1>cells[i].xbegin) & (x1<cells[i].xend) & (y1>cells[i].ybegin) & (y1<cells[i].yend) & (z1>cells[i].zbegin) & (z1<cells[i].zend))
		{
			if(cells[i].p_flag==0)
			{
				cells[i].p_ids = malloc(cells[i].p_ids, 1*sizeof(int)) ;
				cells[i].p_ids[p_flag] = id ;
				p_flag++ ;
				exit(0) ;
			}
			else
			{		
				cells[i].p_ids = realloc(cells[i].p_ids, 1*sizeof(int)) ;
				cells[i].p_ids[p_flag] = id ;
				p_flag++ ;
				exit(0) ;
			}
		}
	}
}

//TODO add nearest neighbor logic (3d modular arithmetic.)
int NearestNeighbor(int index)
{
	Neighbor n = malloc(sizeof(Neighbor));
	n->id = FAILURE_VALUE;
	n->distance = FAILURE_VALUE;
	return n;
}

Neighbor NearestNeighborExhaustive(int index, int myChunkStart, int myChunkEnd)
{
	float distance = 10000000;
	int nearestNeighbor = FAILURE_VALUE;
	for(int i = myChunkStart; i < myChunkEnd; i++)
	{
		if(index == i)
		{
			continue;
		}
		float temp = distance(index, i);
		if(temp < distance)
		{
			distance = temp;
			nearestNeighbor = i;
		}
	}
	Neighbor n = malloc(sizeof(Neighbor));
	n.distance = distance;
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
	if(argc < 6)
	{
		perror("Missing arguments. Exepects 1) inputFile.txt ")
	}
	int numPoints = atoi(argv[1]);
	char* inFileName = argv[2];
	int numQueryPoints = atoi(argv[3]);
	char* queryFileName = argv[4];
	char* outFileName = argv[5];

	

	int mpi_myrank;
    int mpi_commsize;
    
    MPI_Init( &argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

    double startTime;
    double endTime;
    if(mpi_myrank == 0)
    {
    	startTime = MPI_Wtime(); //rank 0 is the timekeeper.
    }
	
	FILE *inFile = fopen(inFileName, "r");
	FILE *queyrFile = fopen(queryFileName, "r");
	FILE *outFile = fopen(outFileName, "w");

    point = malloc(sizeof(point3D) * numPoints);
	for(i = 0; i < numPoints; i++)
  	{
  		float x,y,z;
  		fscanf(inFile, "%lf", x);
  		fscanf(inFile, "%lf", y);
  		fscanf(inFile, "%lf", z);
  		allPoints[i]->x = x;
  		allPoints[i]->y = y;
  		allPoints[i]->z = z;
  	}

  	int chunkSize = (numPoints / mpi_commsize);
  	int myChunkStart = mpi_myrank * chunkSize;

  	//place my points on my grid.
  	for(int i = myChunkStart; i < myChunkStart + chunkSize; i++)
  	{
  		home_cell(i);

	int* queryPoints = malloc(sizeof(int) * numQueryPoints);
	for(i = 0; i < numPoints; i++)
  	{
  		float index;
  		fscanf(inFile, "%d", index);
  		queryPoints[i] = index;
  	}
  	for(int i = 0; i < numQueryPoints; i++) //for every query point (basic version 2 search modes.)
  	{
  		Neighbor n = NearestNeighbor(queryPoints[i]);
  		float* minDistance;
  		*minDistance = n->distance;
  		MPI_AllReduce(minDistance, minDistance, 1, MPI_FLOAT, MPI_MIN);

  		if(distance == FAILURE_VALUE)
  		{
  			Neighbor n = NearestNeighbor(queryPoints[i], myChunkStart, myChunkEnd);
  			*minDistance = n->distance;
  			MPI_AllReduce(minDistance, minDistance, 1, MPI_FLOAT, MPI_MIN);

  			if(*minDistance == n->distance)
  			{
  				//TODO we need some kind of locking for ties.
  				if(distance == FAILURE_VALUE)
  				{
  					if(mpi_myrank == 0)
  					{
  						fputs("%d nearest neighbor: Nonexistent\n",queryPoints[i], outFile);
  					}
  				}
  				else
  				{
  					fputs("%d nearest neighbor: %d\n",queryPoints[i], n->id, outfile)
  				}
  			}
  		}
  	}
	if(mpi_myrank == 0)
    {
    	endTime = MPI_Wtime();
    	double totalTime = endTime-startTime;
    	printf("Total elapsed time was %lf\n", totalTime);
    } 

  	MPI_Finalize();




}

		
