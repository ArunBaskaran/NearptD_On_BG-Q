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
#define FAILURE_VALUE 0.0
float queryX, queryY;

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
char* outFileName ;


//---------------------------------Global Variables---------------------------------------------------//



//----------------------------- Functions -----------------------------------------------------------//

void grid_init()
{
	cell_flag = 0 ; 
	cell_x_limit = 0.0 ;
	cell_y_limit = 0.0 ;
	cell_z_limit = 0.0 ;
	//printf("cells per rank is %d\n", cells_per_rank) ;
	int nx = NX/mpi_commsize ;   
	cell_x_begin = mpi_myrank * nx * dx ;
	cell_x_limit = (mpi_myrank + 1) * nx * dx;
	for(i = 0 ; i < nx ; i++)
	{
		for(j = 0 ; j < NY ; j++)
		{
			for(k = 0 ; k < NZ ; k++)
			{
				//if(cell_flag < cells_per_rank)
				//{
				//printf("rank is %d\n", mpi_myrank) ;
				cells[cell_flag].xbegin = (mpi_myrank*nx*dx) + i*dx ;
				//printf("cells[%d].xbegin = %f\n", cell_flag, cells[cell_flag].xbegin) ;
				cells[cell_flag].xend = (mpi_myrank*nx*dx) + (i+1)*dx ;
				//printf("cells[%d].xend = %f\n", cell_flag, cells[cell_flag].xend) ;
				cells[cell_flag].ybegin = j*dy ;
				cells[cell_flag].yend = (j+1)*dy ;
				cells[cell_flag].zbegin = k*dz ;
				cells[cell_flag].zend = (k+1)*dz ;
				//cells[cell_flag].p_ids = malloc(1*sizeof(struct Cell)) ;
				cells[cell_flag].p_flag = 0 ;
				points_in_cell[cell_flag] = 0 ;
				cell_flag++ ;
				/*if((mpi_myrank+1)*(i+1)*dx > cell_x_limit)
					cell_x_limit = (mpi_myrank+1)*(i+1)*dx ;
				if((j+1)*dy > cell_y_limit)
					cell_y_limit = (j+1)*dy ;
				if((k+1)*dz > cell_x_limit)
					cell_z_limit = (k+1)*dz ; */
				//}
			}
		}
	}
}


/*void readInputFile()
{
	MPI_Offset offset = (int)sizeof(double)*mpi_myrank*numPoints_rank*3 ;
	double x,y,z; 
	for(i = 0; i < numPoints_rank; i++)
  	{
		offset += (int)sizeof(double)*i*3  ;
		MPI_File_read_at_all(fh_input , offset , &x , 1 , MPI_DOUBLE , MPI_STATUS_IGNORE);
		MPI_File_read_at_all(fh_input , offset , &y , 1 , MPI_DOUBLE , MPI_STATUS_IGNORE);
		MPI_File_read_at_all(fh_input , offset , &z , 1 , MPI_DOUBLE , MPI_STATUS_IGNORE);		

  		allPoints[mpi_myrank*numPoints_rank + i].x = x;
  		allPoints[mpi_myrank*numPoints_rank + i].y = y;
  		allPoints[mpi_myrank*numPoints_rank + i].z = z;      

		printf("Read Input files\n") ;
		printf("point oda x is %f\n", allPoints[mpi_myrank*numPoints_rank + i].x) ;
		printf("point oda y is %f\n", allPoints[mpi_myrank*numPoints_rank + i].y) ;
		printf("point oda z is %f\n", allPoints[mpi_myrank*numPoints_rank + i].z) ;
  
  	}
}*/
// , 

//A function that identifies a point's cell, and updates the point's id to the cell's list of points. 
void home_cell(int id)
{ 
	int i ;
	double x1 = allPoints[id].x ;
	double y1 = allPoints[id].y ;
	double z1 = allPoints[id].z ;
	//printf("x1 is %f\n", x1) ;
	//printf("y1 is %f\n", y1) ;
	//printf("z1 is %f\n", z1) ;
	for(int i = 0 ; i < cells_per_rank ; i++)
	{
		//printf("x begin is %f and rank is %d\n", cells[i].xbegin, mpi_myrank) ;
		//printf("x end is %f and rank is %d\n", cells[i].xend, mpi_myrank) ;
		if((x1 >= cells[i].xbegin) & (x1 <= cells[i].xend) & (y1>=cells[i].ybegin) & (y1<cells[i].yend) & (z1>=cells[i].zbegin) & (z1<cells[i].zend))
		{
			//printf("z1 is %f\n", z1) ;
			//printf("Point %d stored in cell %d\n", id, i) ;
			/*if(points_in_cell[i]==0)
			{
				//cells[i].p_ids = malloc(1*sizeof(struct Cell)) ;
				cells[i].p_ids[0] = id ;
				//points_in_cell[i] = 1 ;
				//exit(0) ;
			}*/
			//else
			//{		
				//cells[i].p_ids = realloc(cells[i].p_ids, 1*sizeof(struct Cell)) ;
				cells[i].p_ids[points_in_cell[i]] = id ;
				points_in_cell[i]++ ;
				printf("points in cell %d is %d\n", i,  points_in_cell[i]) ;
				//exit(0) ;
			//}
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
			//printf("Done with get_cell function\n") ;
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
			{printf("point : allPoints[x] = %f and rank is %d\n", allPoints[i].x, mpi_myrank) ; 
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
		FILE *fp;
		fp = fopen("input/input.txt","r");
		fscanf(fp,"%f",&queryX);
		fscanf(fp,"%f",&queryY);
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
	//struct Neighbor n = malloc(sizeof(Neighbor));
	double min_distance = 100000 ;
	//printf("id is %d\n", id) ;
	//printf("In Nearest neighbor function\n") ;
	int success_flag = 0 ;
	int index_new, neighbor_id = 0 ;
	double temp_distance ;
	//printf("Just before querying\n") ;
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
//5: outFile.txt
int main(int argc, char** argv)
{
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
	FILE *inFile = fopen(inFileName, "r");
	FILE *queryFile = fopen(queryFileName, "r");
	FILE *outFile = fopen(outFileName, "w");

        MPI_Init( &argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);
	MPI_Comm file_comm ;
	cells_per_rank = (NX * NY * NZ)/mpi_commsize ;
	cells = malloc(cells_per_rank * sizeof(struct Cell)) ;
	allPoints = malloc(numPoints * sizeof(struct Point_3D)) ;
  	myChunkStart = mpi_myrank * numPoints_rank ;
  	myChunkEnd = (mpi_myrank+1) * numPoints_rank ;
    

        //MPI_File_open(file_comm,inFileName, MPI_MODE_CREATE | MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_input);
//	if(mpi_myrank==0)
//	{
	for(i = 0; i < numPoints; i++)
 	{
   		double x,y,z;
   		fscanf(inFile, "%lf", &x);
   		fscanf(inFile, "%lf", &y);
   		fscanf(inFile, "%lf", &z);
   		allPoints[i].x = x;
   		allPoints[i].y = y;
  		allPoints[i].z = z;
   	}
//	}
	MPI_Barrier(MPI_COMM_WORLD) ;

	grid_init() ;
	//printf("cell_flag is %d\n", cell_flag) ;
	//printf("cells per rank is %d\n", cells_per_rank) ;
	//printf("cell x limit is %f and rank is %d\n", cell_x_limit, mpi_myrank) ;
	assignPointsToCells() ;

	int* queryPoints = malloc(sizeof(int) * numQueryPoints);
	//if(mpi_myrank==0)   // ASSUMING THAT NO. OF RANKS > 1
	//{
	int query_rank = 0 ;
	for(i = 0; i < numQueryPoints; i++)
  	{
  		int index;
  		fscanf(queryFile, "%d", &index);
		if((allPoints[index].x <= cell_x_limit) & (allPoints[index].x > cell_x_begin))
  			{ queryPoints[i] = index;
			  query_rank++ ;
			}
		//printf("query point is 
  	}
	//}
	MPI_Barrier(MPI_COMM_WORLD) ;


	numPoints_rank = numPoints/mpi_commsize ;
	MPI_Barrier(MPI_COMM_WORLD) ;
	MPI_Barrier(MPI_COMM_WORLD) ;

        double startTime;
        double endTime;
        if(mpi_myrank == 0)
        {
    		startTime = MPI_Wtime(); //rank 0 is the timekeeper.
        }


	//printf("Finished assigning\n") ;

	MPI_Barrier(MPI_COMM_WORLD) ;

	int numQueryPoints_rank = numQueryPoints/mpi_commsize ;
	for(int i = 0; i < query_rank ; i++) //for every query point (basic version 2 search modes.)
  	{
   		struct Neighbor n = NearestNeighbor(queryPoints[i]);
   		float minDistance;
		printf("In the nearest neighbor search\n") ;
   		minDistance = n.distance;
		printf("Back from nearest neighbor search\n") ;
		fprintf(outFile , "%d nearest neighbor: %d\n",queryPoints[i], n.id) ;
		//MPI_Barrier(MPI_COMM_WORLD) ;
   		//MPI_Allreduce(&minDistance, &minDistance, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
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


		
