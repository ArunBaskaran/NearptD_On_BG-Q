//Specify file name in main
#include <stdio.h>
#include <string.h>
#include <math.h>

void roundIOvalues(char* fname){
	FILE* fp_input;
	FILE* fp_output;

	fp_input = fopen(fname,"r");
	char outname[15];
	char* extension = "_o";
	strcpy(outname, fname);
	strcpy(outname + strlen(fname), extension);
	fp_output=fopen(outname,"w+");

	float f;
	int count = 0;
	int lines = 0;
	while (!feof(fp_input)){
		fscanf(fp_input, "%f", &f);
		if (count == 2){
					printf("%d\n",lines);

			count = 0;
			lines++;
			fprintf(fp_output, "%.6lf%s", fabs(f), "\n"); 
			continue;
		}
		fprintf(fp_output,"%.6lf%s", fabs(f), " ");
		count ++;
	}


}

int main(){
	roundIOvalues("top3.ply");
}
