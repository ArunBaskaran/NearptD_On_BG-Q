//@authors Evan Maicus, Arun Baskaran, Rahul Divekar, Sidharth Prabhakaran
//The following represents an MPI implementation of the NearptD algorithm.
//See the main function (bottom) for execution instructions
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include <math.h>


#define NX 16
#define NY 16
#define NZ 16
#define dx 0.01
#define dy 0.01
#define dz 0.01
#define FAILURE_VALUE -1
double queryX, queryY, queryZ;

// -----------------------------------Structures-------------------------------------------------------//

//A cell in our world grid. Holds the ids of the points contained within it.
struct Cell
{
	double xbegin, xend, ybegin, yend, zbegin, zend ;
	int * p_ids ;  // A dynamic array to store the ids of points in the cell
	int p_flag ; // An iterative variable that will give the number of points in the cell 
} ;

//A point in 3d space
struct Point_3D
{
	double x, y, z ;
} ;

//A distance and an id. For use when returning nearest neighbor.
struct Neighbor
{
	double distance;
	int id;
};

// -----------------------------------\Structures-------------------------------------------------------//





//---------------------------------Global Variables---------------------------------------------------//

int i, j, k ;
double Lx, Ly, Lz ;
//int cell_flag ;
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
int no_of_cells ;
double cell_x_begin,  cell_x_limit , cell_y_limit, cell_z_limit ;
struct Point_3D * allPoints ;
struct Point_3D * queryPoints ;
struct Cell * cells ;
int cells_per_rank ;

MPI_File fh_input ;
MPI_File fh_query ;
MPI_File fh_output ;

char* inFileName ;
char* queryFileName ;
//char* outFileName ;


//--------------------------------- \Global Variables---------------------------------------------------//



//----------------------------- Functions -----------------------------------------------------------//

//initializes the NX x NY x NZ grid using resolution dx, dy, dz
void grid_init()
{
	int flag = 0 ; 
	for(i = 0 ; i < NX ; i++)
	{
		for(j = 0 ; j < NY ; j++)
		{
			for(k = 0 ; k < NZ ; k++)
			{
				cells[flag].xbegin = i*dx ;
				cells[flag].xend = (i+1)*dx ;
				cells[flag].ybegin = j*dy ;
				cells[flag].yend = (j+1)*dy ;
				cells[flag].zbegin = k*dz ;
				cells[flag].zend = (k+1)*dz ;
				points_in_cell[flag] = 0 ;
			        cells[i].p_ids = malloc(1*sizeof(int)) ;
				flag++ ;

			}
		}
	}
}

//assigns a point to the appropraite cells in grid space.
void home_cell(int id)
{ 
	int i ;

	double x1 = allPoints[id].x ; //get the coordinates of the point
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;

	int xGrid = x1/dx; //scale them into grid space.
	int yGrid = y1/dy;
	int zGrid = z1/dz;

	int index_new = xGrid*NY*NZ + yGrid*NZ + zGrid ; //map them into 1D space.
	//printf("index_new is %d\n", index_new) ;
	cells[index_new].p_ids[points_in_cell[index_new]] = id ; //update
	points_in_cell[index_new]++ ;
	cells[index_new].p_ids = realloc(cells[index_new].p_ids, (points_in_cell[index_new]+1)*sizeof(int)) ;
}


//Given an id, retrieves the proper cell in grid space.
int get_cell(int id)
{
	double x1 = queryPoints[id].x ;
	double y1 = queryPoints[id].y ;
	double z1 = queryPoints[id].z ;
	for(int i = 0 ; i < cell_max ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			return i ;
		}
	}



}

//assigns all points to the appropraite cells in grid space.
void assignPointsToCells()
{
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;

  	for(int i = 0 ; i < numPoints ; i++)
	{
		home_cell(i) ; 		  
	}
}	

//computes the distance between two provided points
double distance(struct Point_3D point1, struct Point_3D point2)
{

	double distance1=0.0;

	distance1 = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2));


	return distance1;
}


/*
Searches cells in immediate vacinity of id to find the nearest neighbor among them. Returns distance 100000 and id -1 if none are found.
Neighboring cells are : index+1, index+1-(NX*NY), index+1+(NX*NY), index+1+NY, index+1-NY, index+1-NY-(NX*NY), index+1-NY+(NX*NY),
 index+1+NY-(NX*NY), index+1+NY-(NX*NY),///// index-1, index-1-(NX*NY), index-1+(NX*NY), index-1+NY, index-1-NY, index-1-NY-(NX*NY), 
 index-1-NY+(NX*NY), index-1+NY-(NX*NY), index-1+NY-(NX*NY), index+NY, index-NY, index + (NX * NY) , index - (NX * NY) 
index-NY-(NX*NY), index-NY+(NX*NY), index+NY-(NX*NY), index+NY-(NX*NY)   (26 cells for 3D)
*/
struct Neighbor NearestNeighbor(int id)
{
	double min_distance = 100000.0 ; //faiure value
	int success_flag = 0 ;
	int index_new, neighbor_id = 0 ;
	double temp_distance ;
	int home_cell_id = get_cell(id) ; //get the starting cell

	//first, we check the immediate cell. This is an optomization to avoid checking adjacent cells if possible.
	index_new = home_cell_id ;
	if((index_new >= 0) & (index_new < cell_max)) //if the index is valid (this should always be true)
 	{ 
 		int temp = 0;
 		while(temp < points_in_cell[index_new]) //while there are more points in the cell
 		{					
 			int temp_id = cells[index_new].p_ids[temp] ;
 			struct Point_3D point1 = queryPoints[id] ;
 			struct Point_3D point2 = allPoints[temp_id] ;
 			temp_distance = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2)) ; //compute the distance
 			if((temp_distance < min_distance) & (temp_distance > 0.0)) //if it is better than our best nearest neighbor (and is not ourself)
 			{
 				min_distance = temp_distance ; //store it
 				neighbor_id = cells[index_new].p_ids[temp] ;
 				success_flag++ ;
 			}
 		
 			temp++ ; 
 		}
 	}
 	if(min_distance != 100000.0) //if we found a neighbor in the immediate cell, return it.
 	{
 		struct Neighbor n ;
 		n.id = neighbor_id ;
 		n.distance = min_distance ;
 		return n;	
 	}

	for(int i = -1 ; i<=1 ; i++) //if we couldn't find a neighbor in our immediate cell, check all adjacent cells.
	{
		for(int j = -1 ; j <= 1 ; j++)
		{
			for(int k = -1 ; k <= 1 ; k++)
			{ 
				index_new = home_cell_id + i + j*NY + k*(NX*NY) ; 
				if(i ==0 && j ==0 && k == 0)
 				{
 					continue;
 				}
				//printf("In index %d\n", index_new) ;
				if((index_new >= 0) & (index_new < cell_max))
				{ 
					int temp = 0;
					while(temp < points_in_cell[index_new])
					{					  
						int temp_id = cells[index_new].p_ids[temp] ;
						struct Point_3D point1 = queryPoints[id] ;
						struct Point_3D point2 = allPoints[temp_id] ;
						temp_distance = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2)) ;
						if((temp_distance < min_distance) & (temp_distance > 0.0))
						{
							min_distance = temp_distance ;
							neighbor_id = cells[index_new].p_ids[temp] ;
							success_flag++ ;
						}
					
						temp++ ; 
					}
				}

			}
		}
	}   // for statement's end brace

	struct Neighbor n ; //return the appropriate neighbor.
	if(success_flag==0)
	{
		n.id = FAILURE_VALUE ;
		n.distance = 100000.0 ;
	}
	else
	{
		n.id = neighbor_id ;
		n.distance = min_distance ;
	}
	return n;	

 




}


//This function does a brute force check of every point on the grid to find the nearest neighbor among them.
struct Neighbor NearestNeighborExhaustive(int index)
{
	double distance_final = 100000.0;
	int nearestNeighbor = FAILURE_VALUE;
	for(int i = 0; i < numPoints_rank; i++)
	{
		if(index == i)
		{
			continue;
		}
		struct Point_3D point1 = queryPoints[index] ;
		struct Point_3D point2 = allPoints[i] ;
		double temp = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2));
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






//--------------------------------MAIN--------------------------------------------------------------//

//args
//1: total # input points
//2: inputFile.txt
//3: total # input query points
//4: queryFile.txt
//5: outFile.txt    //Lets not give outfile 
int main(int argc, char** argv)
{
	if(argc < 4)
	{
		perror("Missing arguments.") ;
	}
	numPoints = atoi(argv[1]);
	inFileName = argv[2];
	numQueryPoints = atoi(argv[3]);
	queryFileName = argv[4];

	char* outFileName = argv[5];

    MPI_Init( &argc, &argv); //initialize mpi.
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);
	Lx = dx * NX ;
	Ly = dy * NY ;
	Lz = dz * NZ ;
 	numPoints_rank = numPoints/mpi_commsize ; //the number of points wer are responsible for.
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;

 	int points_count = 0 ;
	int query_rank = 0 ;
  
    cells = malloc(cell_max * sizeof(struct Cell)) ;
    allPoints = malloc(numPoints_rank * sizeof(struct Point_3D)) ;
    queryPoints = malloc(sizeof(struct Point_3D) * numQueryPoints);

	FILE *inFile = fopen(inFileName, "r");

	fseek(inFile, mpi_myrank*numPoints_rank*27, SEEK_SET) ;

	for(i = 0; i < numPoints_rank; i++) //read in the appropriate points.
 	{

		double x,y,z;
		fscanf(inFile, "%lf", &x);
		fscanf(inFile, "%lf", &y);
		fscanf(inFile, "%lf", &z);
		allPoints[points_count].x = x;
		allPoints[points_count].y = y;
		allPoints[points_count].z = z;
		points_count++ ;
	}

	fclose(inFile) ;
	
	//-----------------READING THE QUERY POINTS----------------//
	FILE *queryFile = fopen(queryFileName, "r");
    for(int i = 0 ; i<numQueryPoints ; i++)
    {
            fscanf(queryFile,"%lf",&queryX);
            fscanf(queryFile,"%lf",&queryY);
            fscanf(queryFile,"%lf",&queryZ);
            queryPoints[i].x = queryX ;
            queryPoints[i].y = queryY ;
            queryPoints[i].z = queryZ ;
	}
	fclose(queryFile) ;
	//-----------------READING THE QUERY POINTS----------------//

    MPI_Barrier(MPI_COMM_WORLD) ;


    double startTime;
    double endTime;
    if(mpi_myrank == 0)
    {
            startTime = MPI_Wtime(); //rank 0 is the timekeeper.
    }


	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//

	//printf("No of points input by rank %d : %d\n", points_count, mpi_myrank);

	for(i = 0 ; i < cell_max ; i++)
	{
	        cells[i].p_ids = malloc(1*sizeof(int)) ;
	}

	grid_init() ;

	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//


	//--------------Assignment of points to cells----------//
	for(i = 0 ; i < points_count ; i++)
	{
		home_cell(i) ;
	}
	//--------------Assignment of points to cells----------//

	if(mpi_myrank == 0)
	{
		FILE *outFile = fopen(outFileName, "w");
		fprintf(outFile, "Beginning run.\n");
		fclose(outFile) ;
	}
	

	int flag = 0 ;
	for(int i = 0; i < numQueryPoints ; i++) //for every query point (basic version 2 search modes.)
  	{
   		struct Neighbor n = NearestNeighbor(i); //see if its nearest neighbor among this ranks points is in the immediate vacinity.
   		double minDistance;
   		minDistance = n.distance;
		double minDistance1 ;
   		MPI_Allreduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); //Reduce to find the minimum nearest neighbor
  		if(minDistance == 100000.0) //if no neighbors were found among all mpi ranks in the immediate vacinity
  		{
			flag++ ;
  			struct Neighbor n1 = NearestNeighborExhaustive(i); //do an exhaustive search of the board
	   		minDistance = n1.distance;
			MPI_Reduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); //and reduce again.
  			if(minDistance1 == 100000)
  			{
				if(mpi_myrank == 0)
				{
					FILE *outFile = fopen(outFileName, "a");
					fprintf(outFile , "Query point no.%d's nearest neighbor: Nonexistent\n", i);
					fclose(outFile) ;
				}
  			}
			else
			{
				if(mpi_myrank == 0)
				{
					FILE *outFile = fopen(outFileName, "a");
					fprintf(outFile , "Query point no.%d's nearest neighbor was found with distance %lf\n",i, minDistance1) ;
					fclose(outFile) ;
				}
			}
  		}
		else
		{
			if(mpi_myrank==0)
			{
				FILE *outFile = fopen(outFileName, "a");
 				fprintf(outFile , "Query point no.%d's nearest neighbor was found at distance %lf\n",i, minDistance1) ;
				fclose(outFile) ;
			}
		}


  	} 

	MPI_Barrier(MPI_COMM_WORLD) ;

	if(mpi_myrank == 0)
    {
		endTime = MPI_Wtime();
		double totalTime = endTime-startTime;
		printf("Total elapsed time was %lf\n", totalTime);
    } 

	free(cells) ;
	free(allPoints) ;


  	MPI_Finalize();

}

//--------------------------------MAIN--------------------------------------------------------------//


		
