#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

unsigned long long SEED = 3001888822891; // have value
unsigned long long RANV;
int RANI = 0;

//declare function
double Ranq1();
void normal(double sig, double *n1, double *n2);
int sgn(double L);
double minabs(double L1, double L2);
double triangle(double L1, double L2);

int main() {
    //now for matrix 7 * 7
    double x, y;               // for normal random variable
    int i, j, k, m;         // for counting
    //int step;
    int num = 0;                // number of block to transmit
    double EbN0 = 0.9844;
    printf("%g\n", EbN0/10);
    double sigma = sqrt(pow(pow(10, EbN0/10.0), -1));
    printf("standard deviation = %g\n",sigma);
    EbN0 = pow(10, EbN0/10.0);
    printf("%g\n", EbN0);

    int n, rc;  // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    
    int *codarray;          //codeword, Dynamic memory allocation, length = codarraylen
    int codarraylen = n;

    //normal(0.9844, &x, &y);
    //printf("%g %g\n", x, y);

    codarray = (int *)malloc(codarraylen * sizeof(int));
    if( codarray == NULL ) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }

    int *codearray;         //0->1; 1->-1,  Dynamic memory allocation, length = codearraylen
    int codearraylen = n;
    
    codearray = (int *)malloc( codearraylen * sizeof(int) );
    if( codearray == NULL ) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }

    double *outp;           // codeword + noise, Dynamic memory allocation, length = outparray
    int outparray = n;

    outp = (double *)malloc( outparray * sizeof(double) );
    if (outp == NULL) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }

    int *output;            // result of interative algorithm decoding, Dynamic memory allocation, length = outputarray
    int outputarray = n;

    output = (int *)malloc( outputarray * sizeof(int) );
    if (output == NULL) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    
    int *L;                 // check node connect 6 variable nodes
    //int Llenrow = 408;
    //int Llencolumn = 6;
    int Llenrow = n;        // n=7
    int Llencolumn = dv;

    L = (int *)malloc(Llenrow * Llencolumn * sizeof(int));
    if (L == NULL) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    
    int *M;
    //int MLlenrow = 816;
    //int Mlencolumn = 3;
    int Mlenrow = rc;       // k=7
    int Mlencolumn = dc;   //dc = 3

    M = (int *)malloc(Mlenrow * Mlencolumn * sizeof(int));
    if (L == NULL) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }

    double *Lj;                // LLR, length = n
    //int n = 7, rc;              // now suppose n = 7, future read data from txt

    Lj = (double *)malloc(n * sizeof(double));
    if (Lj == NULL) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    // DECLARATION ENDING

    // CODING PART NOW
    FILE *fpr;
    fpr=fopen("parity1.txt","r");
    fscanf(fpr,"%d",&n);
    //printf("n = %d\n", n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d\n", n);
    printf("row = %d\n", rc);
    fscanf(fpr,"%d",&dv);
    fscanf(fpr,"%d",&dc);
    printf("dv = %d\n", dv);
    printf("dc = %d\n", dc); 
    for (i = 0; i < n; i++) fscanf(fpr,"%d",&a);
    printf("%d\n", a);
    for (i = 0; i < rc; i++) fscanf(fpr,"%d",&a);
    printf("%d\n", a);
    int  index;
    for (j = 0; j < n; j++) {
        for (i = 0; i < dv; i++) {
            index = j * dv + i;
            fscanf(fpr,"%d",&M[index]);
            printf("M[%d] = %d ", index, M[index]);
        }
        //printf("\n");
    }

    for (i = 0; i < rc; i++) {
        for (j = 0; j < dc; j++) {
            index = i * dc + j;
            fscanf(fpr,"%d",&L[index]);
            printf("L[%d] = %d ", index, L[index]);
        }
        //printf("\n");
    }
    fclose(fpr);

    printf("cod = \n");
    for (i = 0; i < codarraylen; i++) {
        codarray[i] = 0;
        printf("%d ", codarray[i]);
    }
    printf("\n");
    
    printf("code = \n");
    for (i = 0; i < codearraylen; i++) {
        if (codarray[i] == 0) codearray[i] = 1;
        else codearray[i] = -1;
        printf("%d ", codearray[i]);
    }
    printf("\n");

    /*for (i = 0; i < outparray/2; i++) {           // for 816 *408 and 8000 * 4000
        normal(sigma, &x, &y);
        outp[i] = codearray[i] + sigma * x;
        outp[outparray - 1 - i] = codearray[outparray - 1 - i] + sigma * y;
    }*/
    for (i = 0; i < outparray; i++){            // for 7 * 7
        normal(sigma, &x, &y);
        //printf("x[%d] = %g; y[%d] = %g\n", i, x, i, y);
        outp[i] = codearray[i] + sigma * x;
        //printf("outp[%d] =  %g;\n", i, outp[i]);
    }

    printf("outp = \n");
    for(i = 0; i < outparray; i++) {
            printf("%g ", outp[i]);
    }
    printf("\n");

    /*printf("Lj[i] = \n");
    for(i = 0; i < outparray; i++) {
            Lj[i] = 4 * log2(408) * 0.5 * EbN0 * outp[i];
            //Lj[i] =   4 * 0.5 * ebn0 * outp[i];     //  0.5 * 1.2544 = Es/N0
            //printf("%g ", Lj[i]);
        }*/


    free(codarray);
    free(codearray);
    free(outp);
    free(output);
    free(Lj);
    return 0;
}
void normal(double sig, double *n1, double *n2)
{   
    double x1,x2;
    double s;

    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    *n1 = sig * x1 * sqrt((-2.0 * log(s))/ s);
    *n2 = sig * x2 * sqrt((-2.0 * log(s))/ s);
    
}

double Ranq1() {
    if ( RANI == 0 ){
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;

    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

int sgn (double L){
    if (L >= 0) return 1;
    else return -1;
}

double minabs(double L1, double L2) {
    if(L1 <= 0) L1 = (-1) * L1;
    else L1 = L1;
    if(L2 <= 0) L2 = (-1) * L2;
    else L2 = L2;
    if(L1>=L2) return L2;
    else return L1;
}

double triangle(double L1, double L2) {
    double temp1, temp2;
    double ope1,ope2;
    double answer;
    ope1 = L1 + L2;
    ope2 = L1 - L2;
    if (ope1 <= 0) ope1 = ope1;
    else ope1 = (-1) * ope1;
    if (ope2 <= 0) ope2 = ope2;
    else ope2 = (-1) * ope2;
    temp1 = 1 + exp(ope1);
    temp2 = 1 + exp(ope2);
    answer = log(temp1 / temp2);
    return answer;
}