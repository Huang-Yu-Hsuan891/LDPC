#include<stdio.h>
#include<stdlib.h>

int main() {

    int i,j;
    int n,rc;   // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    int input;
    int H[408][816];
    int dv1,dv2,dv3;

    FILE *fpr;
    // open file
    fpr=fopen("parity.txt","r");
    // read file
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d\n", n);
    printf("row = %d\n", rc);

    fscanf(fpr,"%d",&dv);
    fscanf(fpr,"%d",&dc);
    printf("dv = %d\n", dv);
    printf("dc = %d\n", dc);

    for (i = 0; i < 816; i++) fscanf(fpr,"%d",&a);
    printf("dv = %d\n", a);
    for (i = 0; i < 408; i++) fscanf(fpr,"%d",&a);
    printf("dc = %d\n", a);

    // parity check matrix
    for (j = 0; j < 816; j++) {
        fscanf(fpr,"%d",&dv1);
        fscanf(fpr,"%d",&dv2);
        fscanf(fpr,"%d",&dv3);
        //printf("dv1 = %d; dv2 = %d; dv3 = %d",dv1,dv2,dv3);
        for (i = 1; i < 409; i++) {
            if (i == dv1 || i == dv2 || i == dv3) H[i-1][j] = 1;
            else H[i-1][j] = 0;
        }
    }
    for (i = 0; i < 408; i++) {
        for (j = 0; j < 816; j++) {
            printf("%d ", H[i][0]);
        }
        printf("\n");
    }
    fclose(fpr);

    ///////////////////////////
    
    //
    return 0;
}