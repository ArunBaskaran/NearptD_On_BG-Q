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

int i, j, k ;
int Lx, Ly, Lz ;
Lx = dx * NX ;
Ly = dy * NY ;
Lz = dz * NZ ;
int cell_flag = 0 ;
const int cell_max = NX * NY * NZ ;

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

		
