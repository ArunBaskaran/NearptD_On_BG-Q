#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include <math.h>


#define NX 128
#define NY 128
#define NZ 128
#define dx 0.001
#define dy 0.001
#define dz 0.001
#define FAILURE_VALUE 0.0
float queryX, queryY, queryZ;

// -----------------------------------Structures-------------------------------------------------------//

struct Cell
{
	double xbegin, xend, ybegin, yend, zbegin, zend ;
	int * p_ids ;  // A dynamic array to store the ids of points in the cell
	int p_flag ; // An iterative variable that will give the number of points in the cell 
} ;


struct Point_3D
{
	double x, y, z ;
} ;

struct Neighbor
{
	float distance;
	int id;
};

// -----------------------------------Structures-------------------------------------------------------//





//---------------------------------Global Variables---------------------------------------------------//

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
int ** point_ids_cells ;
int points_in_cell[NX * NY * NZ] ;
struct Point_3D * queryPoints ;
int no_of_cells ;
double cell_x_begin,  cell_x_limit , cell_y_limit, cell_z_limit ;
struct Point_3D * allPoints ;
struct Cell * cells ;
int cells_per_rank ;

MPI_File fh_input ;
MPI_File fh_query ;
MPI_File fh_output ;

char* inFileName ;
char* queryFileName ;
//char* outFileName ;


//---------------------------------Global Variables---------------------------------------------------//



//----------------------------- Functions -----------------------------------------------------------//

void grid_init()
{
	cell_flag = 0 ; 
	cell_x_limit = 0.0 ;
	cell_y_limit = 0.0 ;
	cell_z_limit = 0.0 ;
	int nx = NX/mpi_commsize ;   
	cell_x_begin = mpi_myrank * nx * dx ;
	cell_x_limit = (mpi_myrank + 1) * nx * dx;
	for(i = 0 ; i < nx ; i++)
	{
		for(j = 0 ; j < NY ; j++)
		{
			for(k = 0 ; k < NZ ; k++)
			{
				cells[cell_flag].xbegin = (mpi_myrank*nx*dx) + i*dx ;
				cells[cell_flag].xend = (mpi_myrank*nx*dx) + (i+1)*dx ;
				cells[cell_flag].ybegin = j*dy ;
				cells[cell_flag].yend = (j+1)*dy ;
				cells[cell_flag].zbegin = k*dz ;
				cells[cell_flag].zend = (k+1)*dz ;
				cells[cell_flag].p_flag = 0 ;
				points_in_cell[cell_flag] = 0 ;
				cell_flag++ ;

			}
		}
	}
}


void home_cell(int id)
{ 
	int i ;
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;
	for(int i = 0 ; i < cells_per_rank ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
			points_in_cell[i]++ ;
			//printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
}

int get_cell(int id)
{
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;
	for(int i = 0 ; i < cells_per_rank ; i++)
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

  	for(int i = 0 ; i < numPoints ; i++)
  	{
  		if((allPoints[i].x <= cell_x_limit) & (allPoints[i].x > cell_x_begin))    //To DO : Check for points on the borders
			{points_in_cell[i]++ ; }
			  
	}

  	for(int i = 0 ; i < cells_per_rank ; i++)
  	{
  		cells[i].p_ids = malloc(points_in_cell[i]*sizeof(int)) ;
			  
	}

  	for(int i = 0 ; i < numPoints ; i++)
  	{
  		if((allPoints[i].x <= cell_x_limit) & (allPoints[i].x > cell_x_begin))    //To DO : Check for points on the borders
			{//printf("point : allPoints[x] = %f and rank is %d\n", allPoints[i].x, mpi_myrank) ; 
			 home_cell(i) ; }
			  
	}
}	

//Function to read from query input file. Stores in queryX, queryY
//Only rank = 0 reads this
//Input file should be space separated in input/input.txt
void readQueryFile()
{
	if (mpi_myrank == 0)
	{
		for(int i = 0 ; i<numQueryPoints ; i++)
		{
			FILE *fp;
			fp = fopen("input/input.txt","r");
			fscanf(fp,"%f",&queryX);
			fscanf(fp,"%f",&queryY);
			fscanf(fp,"%f",&queryZ);
			queryPoints[i].x = queryX ;
			queryPoints[i].y = queryY ;
			queryPoints[i].z = queryZ ;			
		}
	}
}


double distance(struct Point_3D point1, struct Point_3D point2)
{

	double distance1=0.0;

	distance1 = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2));


	return distance1;
}



struct Neighbor NearestNeighbor(int id)
{
	double min_distance = 100000 ;
	int success_flag = 0 ;
	int index_new, neighbor_id = 0 ;
	double temp_distance ;
	if((allPoints[id].x <= cell_x_limit) & (allPoints[id].x > cell_x_begin))
	{
	//printf("Inside query function\n") ;
	int home_cell_id = get_cell(id) ;
	//printf("Back to Nearest neighbor function\n") ;
	for(int i = -1 ; i<=1 ; i++)
	{
		for(int j = -1 ; j <= 1 ; j++)
		{
			for(int k = -1 ; k <= 1 ; k++)
			{ 
				index_new = home_cell_id + i + j*NY + k*(NX*NY) ; 
				if((index_new >= 0) & (index_new < cell_max))
				{ 
					int temp = 0;
					while(temp <= points_in_cell[index_new])
					{					
						if(points_in_cell[index_new]>0)
						{    
								int temp_id = cells[index_new].p_ids[temp] ;
								struct Point_3D point1 = allPoints[id] ;
								struct Point_3D point2 = allPoints[temp_id] ;
								temp_distance = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2)) ;
								if(temp_distance < min_distance)
								{
									min_distance = temp_distance ;
									neighbor_id = cells[index_new].p_ids[temp] ;
									success_flag++ ;
								}
						}
						temp++ ; 
					}
				}

			}
		}
	}   // for statement's end brace
	}   //if statement's end brace

	struct Neighbor n ;
	if(success_flag==0)
	{
		n.id = 1 ;
		n.distance = 0.0 ;
	}
	else
	{
		n.id = neighbor_id ;
		n.distance = min_distance ;
	}
	return n;	

 


	// Neighboring cells are : index+1, index+1-(NX*NY), index+1+(NX*NY), index+1+NY, index+1-NY, index+1-NY-(NX*NY), index+1-NY+(NX*NY), index+1+NY-(NX*NY), index+1+NY-(NX*NY),///// index-1, index-1-(NX*NY), index-1+(NX*NY), index-1+NY, index-1-NY, index-1-NY-(NX*NY), index-1-NY+(NX*NY), index-1+NY-(NX*NY), index-1+NY-(NX*NY), index+NY, index-NY, index + (NX * NY) , index - (NX * NY) //// index-NY-(NX*NY), index-NY+(NX*NY), index+NY-(NX*NY), index+NY-(NX*NY)   (26 cells for 3D)


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


//----------------------------- Functions -----------------------------------------------------------//




//--------------------------------MAIN--------------------------------------------------------------//

//args
//1: total # input points
//2: inputFile.txt
//3: total # input query points
//4: queryFile.txt
//5: outFile.txt    //Lets not give outfile 
int main(int argc, char** argv)
{
	Lx = dx * NX ;
	Ly = dy * NY ;
	Lz = dz * NZ ;
	cell_flag = 0 ;

	if(argc < 4)
	{
		perror("Missing arguments.") ;
	}
	numPoints = atoi(argv[1]);
	inFileName = argv[2];
	numQueryPoints = atoi(argv[3]);
	queryFileName = argv[4];
	//outFileName = argv[5] ;
	FILE *inFile = fopen(inFileName, "r");
	FILE *queryFile = fopen(queryFileName, "r");
	//FILE *outFile = fopen(outFileName, "w");

	char outFileName[80] = { } ;

        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);
	MPI_Comm file_comm ;
	cells_per_rank = (NX * NY * NZ)/mpi_commsize ;
	cells = malloc(cells_per_rank * sizeof(struct Cell)) ;
	allPoints = malloc(numPoints * sizeof(struct Point_3D)) ;
	queryPoints = malloc(numQueryPoints * sizeof(struct Point_3D)) ;
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;
  

        double startTime;
        double endTime;
        if(mpi_myrank == 0)
        {
    		startTime = MPI_Wtime(); //rank 0 is the timekeeper.
        }
  
	numPoints_rank = numPoints/mpi_commsize ;
	
	fseek(inFile, mpi_myrank*numPoints_rank*30, SEEK_SET) ;

	for(i = 0; i < numPoints_rank; i++)
 	{
   		float x,y,z;
   		fscanf(inFile, "%f", &x);
   		fscanf(inFile, "%f", &y);
   		fscanf(inFile, "%f", &z);
		printf(" x, y, and z are %f, %f, and %f. Rank is %d\n", x, y, z, mpi_myrank) ; 
   		allPoints[i].x = x;
   		allPoints[i].y = y;
  		allPoints[i].z = z;
		//printf("Size of row is %lu\n", sizeof(x)+sizeof(y)+sizeof(z)) ; 
   	}

	MPI_Barrier(MPI_COMM_WORLD) ;

	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//

	grid_init() ;
	assignPointsToCells() ;

	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//


	int* queryPoints = malloc(sizeof(int) * numQueryPoints);
	int query_rank = 0 ;

	//-----------------READING THE QUERY POINTS----------------//

	for(i = 0; i < numQueryPoints; i++)
  	{
  		int index;
  		fscanf(queryFile, "%d", &index);
		if((allPoints[index].x <= cell_x_limit) & (allPoints[index].x > cell_x_begin))
  			{ queryPoints[i] = index;
			  query_rank++ ;
			}
  	}
	MPI_Barrier(MPI_COMM_WORLD) ;
	sprintf(outFileName, "output_%d.txt", mpi_myrank);
	FILE *outFile = fopen(outFileName, "w");

	//-----------------READING THE QUERY POINTS----------------//


	//-----------------NEAREST NEIGHBOR FUNCTION----------------//

	for(int i = 0; i < query_rank ; i++) //for every query point (basic version 2 search modes.)
  	{
   		struct Neighbor n = NearestNeighbor(queryPoints[i]);
   		float minDistance;
		printf("In the nearest neighbor search\n") ;
   		minDistance = n.distance;
		printf("Back from nearest neighbor search\n") ;
		fprintf(outFile , "%d nearest neighbor: %d\n",queryPoints[i], n.id) ;
		MPI_Barrier(MPI_COMM_WORLD) ;
   		MPI_Allreduce(&minDistance, &minDistance, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD) ;
  		if(n.distance == (float)FAILURE_VALUE)
  		{
  			//struct Neighbor n = NearestNeighborExhaustive(queryPoints[i], myChunkStart, myChunkEnd);
  			minDistance = n.distance;
			//MPI_Barrier(MPI_COMM_WORLD) ;
  			//MPI_Allreduce(&minDistance, &minDistance, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

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
  		}


  	}

	//-----------------NEAREST NEIGHBOR FUNCTION----------------//

	MPI_Barrier(MPI_COMM_WORLD) ;

	if(mpi_myrank == 0)
        {
    		endTime = MPI_Wtime();
    		double totalTime = endTime-startTime;
    		printf("Total elapsed time was %lf\n", totalTime);
        } 

	fclose(inFile) ;
	fclose(queryFile) ;
	fclose(outFile) ;

	//MPI_File_close(&fh_input) ;
	//MPI_File_close(&fh_query) ;
	//MPI_File_close(&fh_output) ;
	free(cells) ;
	free(allPoints) ;

  	MPI_Finalize();

}

//--------------------------------MAIN--------------------------------------------------------------//


		
