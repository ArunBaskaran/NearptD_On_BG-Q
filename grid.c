#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include <math.h>



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
double Lx, Ly, Lz ;
int NX, NY, NZ;
//int cell_flag ;
int mpi_myrank;
int mpi_commsize;
int cell_max;
int* points_in_cell;
int numPoints ;
int numPoints_rank ;
int myChunkStart ;
int myChunkEnd ;
int numQueryPoints ;
int ** point_ids_cells ;
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


/*void home_cell(int id)
{ 
	int i ;

	//Reducing the number of iterations that has to be performed
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;
	//printf("Inside home_cell\n") ;
	//printf("x1 is %lf\n", x1) ;
//	Reducing the number of iterations that has to be performed
	if((x1 > 0) && (x1 < Lx/4)) 
	{ printf("Inside first if\n") ;
	for(int i = 0 ; i < cell_max/4 ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

	if((x1 > Lx/4) && (x1 < Lx/2)) 
	{ printf("Inside second if\n") ;
	for(int i = cell_max/4 ; i < 2*cell_max/4 ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

	if((x1 > Lx/2) && (x1 < 3*Lx/4)) 
	{ printf("Inside third if\n") ;
	for(int i = cell_max/2 ; i < 3*cell_max/4 ; i++)
	{
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
										//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

	if((x1 > 3*Lx/4) && (x1 < Lx)) 
	{ printf("Inside fourth if\n") ;
	for(int i = 3*cell_max/4 ; i < cell_max ; i++)
	{	//printf("cell_id is %d\n", i) ;
		if((x1 >= cells[i].xbegin) && (x1 < cells[i].xend) && (y1>=cells[i].ybegin) && (y1<cells[i].yend) && (z1>=cells[i].zbegin) && (z1<cells[i].zend))
		{
			cells[i].p_ids[points_in_cell[i]] = id ;
			printf("Got a match\n") ;							//cells[i].p_ids = realloc(cells[i].p_ids, 1**sizeof(int)) ;
			//printf("In index %d\n", i) ;
			points_in_cell[i]++ ;
			printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
		}
	}
	}

}*/


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
	// printf("%d of %d\n", index_new, cell_max);
	cells[index_new].p_ids[points_in_cell[index_new]] = id ; //update
	points_in_cell[index_new]++ ;

}



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

	//test the i == 0, j== 0, k == 0 case first.

	index_new = home_cell_id ;
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
	if(min_distance != 100000.0)
	{
		struct Neighbor n ;
		n.id = neighbor_id ;
		n.distance = min_distance ;
		return n;	
	}




	for(int i = -1 ; i<=1 ; i++)
	{
		for(int j = -1 ; j <= 1 ; j++)
		{
			for(int k = -1 ; k <= 1 ; k++)
			{ 
				if(i ==0 && j ==0 && k == 0)
				{
					continue;
				}
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
	NX = 1/dx;
	NY = 1/dy;
	NZ = 1/dz;

	cell_max = NX * NY * NZ ;


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
	
	points_in_cell = malloc(sizeof(int) * NX*NY*NZ) ;

	Lx = dx * NX ;
	Ly = dy * NY ;
	Lz = dz * NZ ;
 	numPoints_rank = numPoints/mpi_commsize ;
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;

 	int points_count = 0 ;
	int query_rank = 0 ;
	if(mpi_myrank==0)
	{
		sprintf(outFileName, "output_%d.txt", mpi_myrank);
	}
  
    cells = malloc(cell_max * sizeof(struct Cell)) ;
    allPoints = malloc(numPoints * sizeof(struct Point_3D)) ;
    queryPoints = malloc(sizeof(struct Point_3D) * numQueryPoints);
    printf("%d malloced %lf gigabytes\n", mpi_myrank, (double)((sizeof(struct Point_3D) * numQueryPoints) + (numPoints * sizeof(struct Point_3D)) + (cell_max * sizeof(struct Cell))) / 1073741824);


	fseek(inFile, mpi_myrank*numPoints_rank*27, SEEK_SET) ;

	for(i = 0; i < numPoints_rank; i++)
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

    for(int i = 0 ; i<numQueryPoints ; i++)
    {
	    fscanf(queryFile,"%lf",&queryX);
	    fscanf(queryFile,"%lf",&queryY);
	    fscanf(queryFile,"%lf",&queryZ);
	    queryPoints[i].x = queryX ;
	    queryPoints[i].y = queryY ;
	    queryPoints[i].z = queryZ ;
	}	

    MPI_Barrier(MPI_COMM_WORLD) ;


    double startTime;
    double endTime;
	if(mpi_myrank == 0)
    {
    	endTime = 0;
    	startTime = MPI_Wtime();
    }

	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//


	for(i = 0 ; i < cell_max ; i++)
	{
	        cells[i].p_ids = malloc(500) ;
	}

	grid_init() ;

	

	//-----------------GRID INITIALIZATION AND POINT ALLOCATION----------------//

	//MPI_Barrier(MPI_COMM_WORLD) ;

	//--------------Assignment of points to cells----------//
	for(i = 0 ; i < points_count ; i++)
	{
		//printf("points count is points_count\n") ;
		home_cell(i) ;
	}
	
	//--------------Assignment of points to cells----------//

	//-----------------READING THE QUERY POINTS----------------//

	/*for(int i = 0 ; i<numQueryPoints ; i++)
	{
		fscanf(queryFile,"%lf",&queryX);
		fscanf(queryFile,"%lf",&queryY);
		fscanf(queryFile,"%lf",&queryZ);
		queryPoints[i].x = queryX ;
		queryPoints[i].y = queryY ;
		queryPoints[i].z = queryZ ;
		//printf("querypoint x is %lf\n", queryX) ;			
	}*/

	//MPI_Barrier(MPI_COMM_WORLD) ;

	//-----------------READING THE QUERY POINTS----------------//

	
	int flag = 0 ;
	for(int i = 0; i < numQueryPoints ; i++) //for every query point (basic version 2 search modes.)
  	{
   		struct Neighbor n = NearestNeighbor(i);
   		double minDistance;
   		minDistance = n.distance;
		double minDistance1 ;
   		MPI_Allreduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   		
  		if(minDistance1 == 100000.0)
  		{
  			printf("THIS HAPPENED!\n");
			flag++ ;
  			struct Neighbor n1 = NearestNeighborExhaustive(i);
	   		minDistance = n1.distance;
	   		//MPI_Allreduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			MPI_Reduce(&minDistance, &minDistance1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
			if(mpi_myrank==0)
			{
				FILE *outFile = fopen(outFileName, "a");
  				if(n1.id == FAILURE_VALUE)
  				{
  					if(mpi_myrank == 0)
  					{
  						fprintf(outFile , "Exhuastive Result: Query point no.%d's nearest neighbor: Nonexistent\n", i);
  					}
  				}
  				else
  				{
  					fprintf(outFile , "Exhuastive Result: Query point no.%d's nearest neighbor: %d\n",i, n.id) ;
  				}
  				fclose(outFile) ;
  			}	
			
  			
  		}
		else
		{
			if(mpi_myrank==0)
			{
				FILE *outFile = fopen(outFileName, "a");
 				fprintf(outFile, "QUICK RESULT: Query point no.%d's nearest neighbor: %d\n",i, n.id) ;
				fclose(outFile) ;
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

	fclose(inFile) ;
	fclose(queryFile) ;

	free(cells) ;
	free(allPoints) ;


  	MPI_Finalize();

}

//--------------------------------MAIN--------------------------------------------------------------//


		
