//Specify file name in main
#include <stdio.h>
#include <string.h>

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
	while (!feof(fp_input)){
		fscanf(fp_input, "%f", &f);
		if (count == 2){
			count = 0;
			fprintf(fp_output, "%.8lf%s", f, "\n"); 
			continue;
		}
		fprintf(fp_output,"%.8lf%s", f, " ");
		count ++;
	}


}

int main(){
	roundIOvalues("top3.ply");
}
