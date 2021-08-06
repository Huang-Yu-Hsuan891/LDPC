#include <stdio.h>
#include <math.h>
#include <stdlib.h>

unsigned long long SEED = 389222231;
unsigned long long RANV;
int RANI = 0;

double Ranq1();
void normal(double sig, double *n1, double *n2);

int sgn(double L);
double minabs(double L1, double L2);
double triangle(double L1, double L2);
 
int main() {
    double x,y;
    int i, j, k, m;
    int cod[7];           // codeword
    int code[7];
    double outp[7];       // codeword+noise
    int output[7];     // out of channel
    int num = 0;                // do compute block
    double Lj[7];         // LLR
    double qij[7][7];   // from down to top
    double uij[7][7];   // from top to down
    int L[7][3];          // check L(i)
    int M[7][3];          // check M(j)
    int n,rc;   // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    int input;
    int H[7][7];
    int Hcheck[7][7];
    double tempqij[2]; 
    double tempuij;
    double temp1uij[3];
    double temp1qij;    
    double qj[7];
    int checkbit[7];
    int s = 0;      // receive 100 error block
    int restart = 0;
    int totalerror=0; 
    int error;
    double sigma;
    double ebn0;

    FILE *fpr;
    fpr=fopen("parity77.txt","r");

    // parity check matrix
    for (j = 0; j < 7; j++) {
        fscanf(fpr,"%d",&M[j][0]);
        fscanf(fpr,"%d",&M[j][1]);
        fscanf(fpr,"%d",&M[j][2]);
        for (i = 1; i < 8; i++) {
            if (i == M[j][0] || i == M[j][1] || i == M[j][2]) H[i-1][j] = 1;
            else H[i-1][j] = 0;
        }
    }
    // check parity check maatrix is same
    for (i = 0; i < 7; i++) {
        fscanf(fpr,"%d",&L[i][0]);
        fscanf(fpr,"%d",&L[i][1]);
        fscanf(fpr,"%d",&L[i][2]);
        for (j = 1; j < 8; j++) {
            if (j == L[i][0] || j == L[i][1] || j == L[i][2] || j == L[i][3] || j == L[i][4] || j == L[i][5]) Hcheck[i][j-1] = 1;
            else Hcheck[i][j-1] = 0;
        }
    }
    fclose(fpr);


    while (s != 98/*s < 100*/) {
        num++;   //compute the number of transmit block 
        // pretend encoder
        printf("\n");
        printf("cod = ");
        for(i = 0; i < 7; i++) {
            cod[i] = 0;         // message
            printf("%d ", cod[i]);
        }
        printf("\n");

        //input to AWGN channel normalized to +-1
        printf("code = ");
        for(i = 0; i < 7; i++) {
            if(cod[i] == 0) code[i] = 1;
            else code[i] = -1;
            printf("%d", code[i]);
        }
        printf("\n");

        ebn0 = 0.9844;
        printf("pow = %g\n",pow(10, ebn0/10));
        
        //s++; 
        
        sigma = sqrt(1.0 / pow(10, ebn0/10));
        printf("sigma = %g\n", sigma);

        // add a gaussian random varible of mean = 0 and variance = sigma ^ 2
        for(i = 0; i < 7; i++) {
            normal(sigma, &x, &y);
            //printf("x[%d] = %g; y[%d] = %g\n", i, x, i, y);
            outp[i] = code[i] + x;
            //outp[815-i] = code[815-i] + y;
            //printf("codeword =outp[%d] =  %g\n", i, outp[i]);
        }
        for(i = 0; i < 7; i++) {
            printf("outp[%d] = %g ", i, outp[i]);
        }
        printf("\n");

        ebn0 = pow(10, ebn0/10);
        printf("ebn0 = %g\n", ebn0);

        //printf("Lj[i] = \n");
        for(i = 0; i < 7; i++) {
            Lj[i] = 4 * 0.5 * ebn0 * outp[i];     //  0.5 * 1.2544 = Es/N0
            //printf("%g ", Lj[i]);
        }
        //printf("\n");
        
        // the interative decoding algotrithm
        // initialization
        //printf("initial qij = ");
        for (j = 0; j < 7; j++) {
            for (i = 0; i < 7; i++) {
                if (H[i][j] == 1) {
                    qij[i][j] = Lj[j];
                    //printf("qij[%d][%d] = %g ", i, j, qij[i][j]);
                }
            }   
        }
        //printf("\n");

        // message passing
        for (k = 0; k < 100/*100*/ && restart != 7; k++) {         // for predetermined
            restart = 0;
            // bottom-up
            for (i = 0; i < 2; i++) {
                tempqij[i] = 0.0;
                //printf("tempqij[%d] = %g\n", i, tempqij[i]);
            }
            //printf("\n");
            int valL;
            int valL2;
            for (i = 0; i < 7; i++) {
                for (j = 0; j < 3; j++) {
                    for (m = 0; m < 2; m++) {
                        if (m < j) {
                            valL = L[i][m]-1;
                            //printf("valL = %d\n", valL);
                            tempqij[m] = qij[i][valL];
                            //printf("qij[%d][%d] = %g; tempqij[%d] = %g\n", i, valL, qij[i][valL], m, tempqij[m]);
                        } 
                        else if (m >= j) {
                            valL = L[i][m+1]-1;
                            //printf("valL = %d\n", valL);
                            tempqij[m] = qij[i][valL];
                            //printf("qij[%d][%d] = %g; tempqij[%d] = %g\n", i, valL, qij[i][valL], m, tempqij[m]);
                        }
                    }
                    tempuij = tempqij[0];
                    for(m = 1; m < 2; m++) {
                        tempuij = sgn(tempuij) * sgn(tempqij[m]) * minabs(tempuij,tempqij[m]) + triangle(tempuij,tempqij[m]); 
                        //printf("tempuij = %g\n", tempuij);
                    }
                    valL2 = L[i][j]-1;
                    uij[i][valL2] = tempuij;
                    //printf("uij[%d][%d] = %g\n",i, valL2, uij[i][valL2]);
                }
            }
            //printf("\n\n");

            // top-down
            for(i = 0; i < 3; i++) {
                temp1uij[i] = 0.0;
                //printf("temp1uij[%d] = %g\n", i, temp1uij[i]);
            }
            for (j = 0; j < 7; j++) {
                for (i = 0; i < 3; i++) {
                    for (m = 0; m < 2; m++) {
                        if (m < i) { 
                            valL = M[j][m] - 1;
                            //printf("valL = %d\n", valL);
                            temp1uij[m] = uij[valL][j];
                            //printf("uij[%d][%d] = %g; temp1uij[%d] = %g\n", valL, j, uij[valL][j], m, temp1uij[m]); 
                        }
                        else if (m >= j) {
                            valL = M[j][m+1] - 1;
                            //printf("valL = %d\n", valL);
                            temp1uij[m] = uij[valL][j];
                            //printf("uij[%d][%d] = %g; temp1uij[%d] = %g\n", valL, j, uij[valL][j], m, temp1uij[m]);
                        }
                    }
                    temp1uij[2] = Lj[j];
                    //printf("temp1uij[2] = %g\n",temp1uij[2]);
                    valL2 = M[j][i] - 1;
                    qij[valL2][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                    //printf("qij[%d][%d] = %g\n",valL2, j, qij[valL2][j]);
                }
            }
            //printf("\n\n");

            // decision
            //printf("output = ");
            for (j = 0; j < 7; j++) {
                qj[j] = Lj[j] + uij[M[j][0]-1][j] + uij[M[j][1]-1][j] + uij[M[j][2]-1][j]; 
                //printf("qj[%d] = %g; Lj[%d] = %g; uij[%d][%d] = %g\n",j, qj[j], j, Lj[j], M[j][0]-1, j, uij[M[j][0]-1][j]);
                //printf("uij[%d][%d] = %g; uij[%d][%d] = %g\n", M[j][1]-1, j, uij[M[j][1]-1][j], M[j][2]-1, j, uij[M[j][2]-1][j]);
                if (qj[j] >= 0) output[j] = 0;
                else if (qj[j] < 0) output[j] = 1;
                //printf("output[%d] = %d \n", j, output[j]);
            }
            //printf("\n\n");

            // to check Hx=0
            //printf("checkbit = \n");            
            for (i = 0; i < 7; i++) {
                checkbit[i] = 0;
                for (j = 0; j < 7; j++) {
                    checkbit[i] += (H[i][j] * output[j]);
                    //printf("checkbit[%d] = %d\n", i, checkbit[i]);
                }
                //printf("checkbit[%d] = %d\n", i, checkbit[i]);
                checkbit[i] = checkbit[i] % 2;
                //printf("checkbit/2[%d] = %d\n", i, checkbit[i]);
            }
            //printf("\n");

            for (i = 0; i < 7; i++) {
                if (checkbit[i] == 0) {
                    restart += 1; // restart = 408 is success
                    //printf("restart[%d] = %d\n", i, restart);
                }
            }
            //printf("restart = %d\n", restart);
            //printf("\n\n");
            if (k == 99 && restart != 408) {
                printf("failure\n");
                s++;
            }
            //restart = 0;
        }
        printf("k[%d] = %d\n", num, k);
        error = 0;
        for(i = 0; i < 7; i++) {
            if (output[i] != cod[i]) {
                error += 1;
                printf("ouput[%d] = %d\n",i, output[i]);
            }
        }
        if (error != 0) s++;
        restart = 0;
        printf("error = %d\n", error);       
        totalerror += error;
    }
    printf("totalerror = %d\n", totalerror);
    printf("num = %d\n", num);
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