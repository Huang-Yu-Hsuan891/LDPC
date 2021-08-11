#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
int main() {
    int i, j, k, m;         // for counting
    int M[816][3];          // check M(j)
    int n,rc;   // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    int input;
    int H[408][816];

    FILE *fpr;

    fpr=fopen("parity.txt","r");
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d\n", n);
    printf("row = %d\n", rc);
    fscanf(fpr,"%d",&dv);
    fscanf(fpr,"%d",&dc);
    printf("dv = %d\n", dv);
    printf("dc = %d\n", dc);

    for (i = 0; i < 816; i++) fscanf(fpr,"%d",&a);
    for (i = 0; i < 408; i++) fscanf(fpr,"%d",&a);

    // parity check matrix
    for (j = 0; j < 816; j++) {
        fscanf(fpr,"%d",&M[j][0]);
        fscanf(fpr,"%d",&M[j][1]);
        fscanf(fpr,"%d",&M[j][2]);
        for (i = 1; i < 409; i++) {
            if (i == M[j][0] || i == M[j][1] || i == M[j][2]) H[i-1][j] = 1;
            else H[i-1][j] = 0;
        }
    }
    // check parity check maatrix is same
    fclose(fpr);
    FILE *outfp;
    
    outfp = fopen("parityout.txt","w");
    for (i = 0; i < 408; i++) {
        for (j = 0; j < 816; j++) {
            fprintf(outfp,"%d ",H[i][j]);
        }
        fprintf(outfp,"\n");
    }
    fclose(outfp);
    return 0;
}