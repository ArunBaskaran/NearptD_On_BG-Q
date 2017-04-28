#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include <math.h>


#define NX 12
#define NY 12
#define NZ 12
#define dx 0.01
#define dy 0.01
#define dz 0.01
#define FAILURE_VALUE -1
double queryX, queryY, queryZ;

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
	double distance;
	int id;
};

// -----------------------------------Structures-------------------------------------------------------//





//---------------------------------Global Variables---------------------------------------------------//

int i, j, k ;
int Lx, Ly, Lz ;
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


//---------------------------------Global Variables---------------------------------------------------//



//----------------------------- Functions -----------------------------------------------------------//

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
				//cells[cell_flag].p_flag = 0 ;
				points_in_cell[flag] = 0 ;
				flag++ ;

			}
		}
	}
}


void home_cell(int id)
{ 
	int i ;

	//Reducing the number of iterations that has to be performed
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;

//	Reducing the number of iterations that has to be performed
	if((x1 > 0) & (x1 < Lx/4)) 
	{
	for(int i = 0 ; i < cell_max/4 ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			//printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

	if((x1 > Lx/4) & (x1 < Lx/2)) 
	{
	for(int i = cell_max/4 ; i < 2*cell_max/4 ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			//printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

	if((x1 > Lx/2) & (x1 < 3*Lx/4)) 
	{
	for(int i = cell_max/2 ; i < 3*cell_max/4 ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			//printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

	if((x1 > 3*Lx/4) & (x1 < Lx)) 
	{
	for(int i = 3*cell_max/4 ; i < cell_max ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			//printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

}


/*void home_cell(int id)
{ 
	int i ;

	//Reducing the number of iterations that has to be performed
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;

//	Reducing the number of iterations that has to be performed

	for(int i = 0 ; i < cell_max ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			//printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}

}*/



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


void assignPointsToCells()
{
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;

  	for(int i = 0 ; i < numPoints ; i++)
	{
			 home_cell(i) ; 
			  
	}
}	

//Function to read from query input file. Stores in queryX, queryY
//Only rank = 0 reads this
//Input file should be space separated in input/input.txt


double distance(struct Point_3D point1, struct Point_3D point2)
{

	double distance1=0.0;

	distance1 = sqrt(pow((point2.x-point1.x),2) + pow((point2.y-point1.y),2) + pow((point2.z-point1.z),2));


	return distance1;
}



struct Neighbor NearestNeighbor(int id)
{
	double min_distance = 100000.0 ;
	int success_flag = 0 ;
	int index_new, neighbor_id = 0 ;
	double temp_distance ;
	int home_cell_id = get_cell(id) ;
	//printf("Cell id is %d\n", home_cell_id) ;
	//printf("Number of points in cell %d is %d\n", home_cell_id, points_in_cell[home_cell_id]) ;
	for(int i = -1 ; i<=1 ; i++)
	{
		for(int j = -1 ; j <= 1 ; j++)
		{
			for(int k = -1 ; k <= 1 ; k++)
			{ 
				index_new = home_cell_id + i + j*NY + k*(NX*NY) ; 
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
						if(temp_distance < min_distance)
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

	struct Neighbor n ;
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

 


	// Neighboring cells are : index+1, index+1-(NX*NY), index+1+(NX*NY), index+1+NY, index+1-NY, index+1-NY-(NX*NY), index+1-NY+(NX*NY), index+1+NY-(NX*NY), index+1+NY-(NX*NY),///// index-1, index-1-(NX*NY), index-1+(NX*NY), index-1+NY, index-1-NY, index-1-NY-(NX*NY), index-1-NY+(NX*NY), index-1+NY-(NX*NY), index-1+NY-(NX*NY), index+NY, index-NY, index + (NX * NY) , index - (NX * NY) //// index-NY-(NX*NY), index-NY+(NX*NY), index+NY-(NX*NY), index+NY-(NX*NY)   (26 cells for 3D)


}


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
	if(argc < 4)
	{
		perror("Missing arguments.") ;
	}
	numPoints = atoi(argv[1]);
	inFileName = argv[2];
	numQueryPoints = atoi(argv[3]);
	queryFileName = argv[4];
	FILE *inFile = fopen(inFileName, "r");
	FILE *queryFile = fopen(queryFileName, "r");

	char outFileName[80] = { } ;

        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

	Lx = dx * NX ;
	Ly = dy * NY ;
	Lz = dz * NZ ;
 	numPoints_rank = numPoints/mpi_commsize ;
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;

 	int points_count = 0 ;
	int query_rank = 0 ;
	sprintf(outFileName, "output_%d.txt", mpi_myrank);
	FILE *outFile = fopen(outFileName, "w");
  
 /*       double startTime;
        double endTime;
        if(mpi_myrank == 0)
        {
    		startTime = MPI_Wtime(); //rank 0 is the timekeeper.
        }
*/
        cells = malloc(cell_max * sizeof(struct Cell)) ;
        allPoints = malloc(numPoints * sizeof(struct Point_3D)) ;
        queryPoints = malloc(sizeof(struct Point_3D) * numQueryPoints);

	fseek(inFile, mpi_myrank*numPoints_rank*27, SEEK_SET) ;

	for(i = 0; i < numPoints_rank; i++)
 	{
		//if((i >= mpi_myrank*numPoints_rank) & (i < (mpi_myrank+1)*numPoints_rank))
	//{
   			double x,y,z;
   			fscanf(inFile, "%lf", &x);
	   		fscanf(inFile, "%lf", &y);
   			fscanf(inFile, "%lf", &z);
			//printf("x , y, z are %lf, %lf, %lf in rank %d\n", x,y,z,mpi_myrank) ;
   			allPoints[points_count].x = x;
   			allPoints[points_count].y = y;
  			allPoints[points_count].z = z;
			//home_cell(points_count) ;
			points_count++ ;
  		//}
	}

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
	        cells[i].p_ids = malloc(points_count*sizeof(int)) ;
	}

	grid_init() ;

	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//

	//MPI_Barrier(MPI_COMM_WORLD) ;

	//--------------Assignment of points to cells----------//
	for(i = 0 ; i < points_count ; i++)
	{
		home_cell(i) ;
	}
	//--------------Assignment of points to cells----------//

	//-----------------READING THE QUERY POINTS----------------//

	for(int i = 0 ; i<numQueryPoints ; i++)
	{
		fscanf(queryFile,"%lf",&queryX);
		fscanf(queryFile,"%lf",&queryY);
		fscanf(queryFile,"%lf",&queryZ);
		queryPoints[i].x = queryX ;
		queryPoints[i].y = queryY ;
		queryPoints[i].z = queryZ ;
		//printf("querypoint x is %lf\n", queryX) ;			
	}

	//MPI_Barrier(MPI_COMM_WORLD) ;

	//-----------------READING THE QUERY POINTS----------------//
	int flag = 0 ;
	for(int i = 0; i < numQueryPoints ; i++) //for every query point (basic version 2 search modes.)
  	{
   		struct Neighbor n = NearestNeighbor(i);
   		double minDistance;
		//printf("In the nearest neighbor search\n") ;
   		minDistance = n.distance;
		printf("minDistance for rank %d : %lf\n", mpi_myrank, n.distance);
		double minDistance1 ;
		//printf("Back from nearest neighbor search\n") ;
		fprintf(outFile , "Query point no.%d's nearest neighbor : %d and the distance is %lf\n",i, n.id, n.distance) ;
	//	MPI_Barrier(MPI_COMM_WORLD) ;
   		MPI_Allreduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD) ;
		//int flag =0 ;
  		if(minDistance1 == 100000.0)
  		{
			flag++ ;
  			struct Neighbor n1 = NearestNeighborExhaustive(i);
	   		minDistance = n1.distance;
	//		MPI_Barrier(MPI_COMM_WORLD) ;
	   		MPI_Allreduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  			//minDistance = n.distance;
  			if(minDistance1 == n1.distance)
  			{
  				//TODO we need some kind of locking for ties.
  				if(n1.id == FAILURE_VALUE)
  				{
  					if(mpi_myrank == 0)
  					{
  						fprintf(outFile , "Query point no.%d's nearest neighbor: Nonexistent\n", i);
  					}
  				}
  				else
  				{
  					fprintf(outFile , "Query point no.%d's nearest neighbor: %d\n",i, n.id) ;
  				}
  			}
  		}


  	} 

	printf("No. of exhaustive searches done by rank %d : %d \n", mpi_myrank, flag) ;

	MPI_Barrier(MPI_COMM_WORLD) ;

	if(mpi_myrank == 0)
        {
    		endTime = MPI_Wtime();
    		double totalTime = endTime-startTime;
    		printf("Total elapsed time was %lf\n", totalTime);
        } 

//	fclose(inFile) ;
//	fclose(queryFile) ;
//	fclose(outFile) ;

	free(cells) ;
	free(allPoints) ;


  	MPI_Finalize();

}

//--------------------------------MAIN--------------------------------------------------------------//


		
