#include <stdio.h>
#include <math.h>
#include <stdlib.h>

unsigned long long SEED = 388891;
unsigned long long RANV;
int RANI = 0;

double Ranq1();
void normal (double sig, double *n1, double *n2);

int sgn(int L);
double minabs(double L1, double L2);
double triangle(double L1, double L2);

int main() {
    double x,y;
    int i, j, k, m;
    int cod[816];           // codeword
    int code[816];
    double outp[816];       // codeword+noise
    int output[816];     // out of channel
    int num;                // do 10000 block

    double Lj[816];         // LLR
    double qij[408][816];   // from down to top
    double uij[408][816];   // from top to down
    int L[408][6];          // check L(i)
    int M[816][3];          // check M(j)

    int n,rc;   // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    int input;
    int H[408][816];
    int Hcheck[408][816];
    double tempqij[5]; 
    double tempuij;

    double temp1uij[3];
    double temp1qij;
    
    double qj[816];
    int checkbit[408];

    int s = 0;
    int restart = 0;
    int totalerror=0; 
    int error;
    double sigma;
    double ebn0;

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
        fscanf(fpr,"%d",&M[j][0]);
        fscanf(fpr,"%d",&M[j][1]);
        fscanf(fpr,"%d",&M[j][2]);
        //printf("dv1 = %d; dv2 = %d; dv3 = %d",dv1,dv2,dv3);
        for (i = 1; i < 409; i++) {
            if (i == M[j][0] || i == M[j][1] || i == M[j][2]) H[i-1][j] = 1;
            else H[i-1][j] = 0;
        }
    }
    // check parity check maatrix is same
    for (i = 0; i < 408; i++) {
        fscanf(fpr,"%d",&L[i][0]);
        fscanf(fpr,"%d",&L[i][1]);
        fscanf(fpr,"%d",&L[i][2]);
        fscanf(fpr,"%d",&L[i][3]);
        fscanf(fpr,"%d",&L[i][4]);
        fscanf(fpr,"%d",&L[i][5]);
        for (j = 1; j < 817; j++) {
            if (j == L[i][0] || j == L[i][1] || j == L[i][2] || j == L[i][3] || j == L[i][4] || j == L[i][5]) Hcheck[i][j-1] = 1;
            else Hcheck[i][j-1] = 0;
        }
    }
    fclose(fpr);


    for (num = 0; num < 10/*10000*/; num++) {
        // pretend encoder
        //SEED = SEED + 100;
        //printf("SEED = %lld", SEED);
        printf("\n");
        for(i = 0; i < 816; i++) {
            cod[i] = 0;         // message
        }

        //input to AWGN channel normalized to +-1
        for(i = 0; i < 816; i++) {
            if(cod[i] == 0) code[i] = 1;
            else code[i] = -1;
        }

        ebn0 = 0.9844;
        sigma = sqrt(1.0 / pow(10, ebn0/10));
        // add a gaussian random varible of mean = 0 and variance = sigma ^ 2
        for(i = 0; i < 408; i++) {
            normal(sigma, &x, &y);
            outp[i] = code[i] + x;
            outp[815-i] = code[815-i] + y;
            //printf("codeword = %g %g \n", outp[i], outp[815-i]);
        }

        //printf("Lj[i] = \n");
        ebn0 = pow(10, ebn0/10);
        for(i = 0; i < 816; i++) {
            Lj[i] = 4 * 0.5 * ebn0 * outp[i];     //  0.5 * 1.2544 = Es/N0
            //printf("%g ", Lj[i]);
        }
        //printf("\n");
        
        // the interative decoding algotrithm
        // initialization
        for (j = 0; j < 816; j++) {
            for (i = 0; i < 408; i++) {
                if (H[i][j] == 1) {
                    qij[i][j] = Lj[j];
                    //printf("qij[%d][%d] = %g", i, j, qij[i][j]);
                }
            }   
        }

        // message passing
        for (k = 0; k < 100 && restart != 408; k++) {         // for predetermined
            restart = 0;
            // bottom-up
            for (i = 0; i < 408; i++) {
                for (j = 0; j < 6; j++) {
                    for (m = 0; m < 5; m++) {
                        if (m < j) {
                        //tempqij[m] = qij[i][M[i][m]-1];
                            tempqij[m] = qij[i][L[i][m]-1];
                        } 
                        else if (m >= j) {
                            //tempqij[m] = qij[i][M[i][m + 1]-1];
                            tempqij[m] = qij[i][L[i][m + 1]-1];
                        }
                    }
                    tempuij = tempqij[0];
                    for(m = 1; m < 5; m++) {
                        tempuij = sgn(tempuij) * sgn(tempqij[m]) * minabs(tempuij,tempqij[m]) + triangle(tempuij,tempqij[m]); 
                    }
                    //uij[i][M[i][j]-1] = tempuij;
                    uij[i][L[i][j]-1] = tempuij;
                }
            }

            // top-down
            for (j = 0; j < 816; j++) {
                for (i = 0; i < 3; i++) {
                    for (m = 0; m < 2; m++) {
                        if (m < i) {
                        //temp1uij[m] = uij[L[m][j] - 1][j]; 
                            temp1uij[m] = uij[M[j][m] - 1][j]; 
                        }
                        else if (m >= j) {
                            //temp1uij[m] = uij[L[m+1][j] - 1][j];
                            temp1uij[m] = uij[M[j][m+1] - 1][j];
                        }
                    }
                    temp1uij[2] = Lj[j];
                    //q[L[i][j] - 1][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                    qij[M[j][i] - 1][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                }
            }

            // decision
            printf("output = ");
            for (j = 0; j < 816; j++) {
                qj[j] = Lj[j] + uij[M[j][0]-1][j] + uij[M[j][1]-1][j] + uij[M[j][2]-1][j]; 
                if (qj[j] >= 0) output[j] = 0;
                else if (qj[j] < 0) output[j] = 1;
                printf("%d ", output[j]);
            }
            printf("\n");

            // to check Hx=0
            printf("checkbit = ");            
            for (i = 0; i < 408; i++) {
                checkbit[i] = 0;
                for (j = 0; j < 816; j++) {
                    checkbit[i] += (H[i][j] * output[j]);
                }
                checkbit[i] = checkbit[i] % 2;
                printf("%d",checkbit[i]);
            }
            printf("\n");

            for (i = 0; i < 408; i++) {
                if (checkbit[i] == 0) restart += 1; // restart = 408 is success
            }
            if (k == 99 && restart != 408) printf("failure\n");
        }
        printf("k[%d] = %d\n", num, k);
        error = 0;
        for(i = 0; i < 816; i++) {
            if (output[i] != cod[i]) {
                error += 1;
                printf("ouput[%d] = %d",i, output[i]);
            }
        }
        restart = 0;
        printf("error = %d\n", error);       
        totalerror += error;
    }
    printf("totalerror = %d\n", totalerror);
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
    *n1 = sig * x1 * sqrt(-2 * log(s) / s);
    *n2 = sig * x2 * sqrt(-2 * log(s) / s);
    
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

int sgn (int L){
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