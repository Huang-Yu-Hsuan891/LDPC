#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


unsigned long long SEED = 3122891; // have value
//unsigned long long SEED = 3881111387891;
unsigned long long RANV;
int RANI = 0;

double Ranq1();
void normal(double sig, double *n1, double *n2);

int sgn(double L);
double minabs(double L1, double L2);
double triangle(double L1, double L2);

int main(){
    double x,y;             // for normal random variable
    int i, j, k, m;         // for counting
    int step;
    int cod[816];           // codeword
    int code[816];
    double outp[816];       // codeword+noise
    int output[816];     // out of channel
    int num = 0;                // do compute block
    double Lj[816];         // LLR
    double qij[3][816];   // from down to top
    //double uij[3][8000];
    double uij[408][6];   // from top to down
    
    int L[408][6];          // check L(i)
    int M[816][3];          // check M(j)
    int n,rc;   // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    int input;
    //int H[4000][8000];
    //int Hcheck[4000][8000];
    double tempqij[5]; 
    double tempuij;
    double temp1uij[3];
    double temp1qij;   
    double qj[816];
    int checkbit[408];
    int s = 0;      // receive 100 error block
    int restart = 0;
    int totalerror=0; 
    int error;
    double sigma;
    double ebn0;
    double ebn0s[6];
    double bers[6];
    double berscompare[6];
    int stp;

    ebn0s[0] = 0.9844;
    ebn0s[1] = 1.289;
    ebn0s[2] = 1.584;
    ebn0s[3] = 1.868;
    ebn0s[4] = 2.411;
    ebn0s[5] = 3.046;

    berscompare[0] = 0.0583;
    berscompare[1] = 0.03052;
    berscompare[2] = 0.01011;
    berscompare[3] = 0.003;
    berscompare[4] = 8.34 * pow(10,-5);
    berscompare[5] = 2.799 * pow(10,-7);

    FILE *fpr;
    fpr=fopen("parity1.txt","r");
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

    for (j = 0; j < 816; j++) {
        fscanf(fpr,"%d",&M[j][0]);
        fscanf(fpr,"%d",&M[j][1]);
        fscanf(fpr,"%d",&M[j][2]);
    }

    for (i = 0; i < 408; i++) {
        fscanf(fpr,"%d",&L[i][0]);
        fscanf(fpr,"%d",&L[i][1]);
        fscanf(fpr,"%d",&L[i][2]);
        fscanf(fpr,"%d",&L[i][3]);
        fscanf(fpr,"%d",&L[i][4]);
        fscanf(fpr,"%d",&L[i][5]);
    }
    fclose(fpr);    

    for(step = 0; step < 6; step++) {
        s = 0;
        num = 0;
        totalerror = 0; 

        while (s < 100/*s < 100*/) {
            num++;   //compute the number of transmit block 
            // pretend encoder
            //printf("\n");
            //printf("cod = ");
            for(i = 0; i < 816; i++) {
                cod[i] = 0;         // message
                //printf("%d ", cod[i]);
            }
            //printf("\n");

            //input to AWGN channel normalized to +-1
            //printf("code = ");
            for(i = 0; i < 816; i++) {
                if(cod[i] == 0) code[i] = 1;
                else code[i] = -1;
                //printf("%d", code[i]);
            }
            //printf("\n");

            ebn0 = ebn0s[step];
            //printf("ebn0s[%d] = %g\n", step, ebn0);
            //printf("pow = %g\n",pow(10, ebn0/10)); 
            sigma = sqrt(1.0 / pow(10, ebn0/10));
            //printf("sigma = %g\n", sigma);

            // add a gaussian random varible of mean = 0 and variance = sigma ^ 2
            for(i = 0; i < 408; i++) {
                normal(sigma, &x, &y);
                //printf("x[%d] = %g; y[%d] = %g\n", i, x, i, y);
                outp[i] = code[i] + /*sigma **/ x;
                outp[815-i] = code[815-i] + /*sigma **/ y;
                //printf("codeword =outp[%d] =  %g;outp[%d] = %g \n", i, outp[i], 815-i, outp[7999-i]);
            }

            //printf("Lj[i] = \n");
            ebn0 = pow(10, ebn0/10);
            for(i = 0; i < 816; i++) {
                Lj[i] = /*log2(408) **/ 4 * 0.5 * ebn0 * outp[i];     //  0.5 * 1.2544 = Es/N0
                //printf("%g ", Lj[i]);
            }
            //printf("\n");
        
            // the interative decoding algotrithm
            // initialization
            for (j = 0; j < 816; j++) {
                for (i = 0; i < 3; i++) {
                        qij[i][j] = Lj[j];
                        //printf("qij[%d][%d] = %g  ", i, j, qij[i][j]);
                }   
            }
            printf("\n");

            // message passing
            for (k = 0; k < 100 && restart != 408; k++) {         // for predetermined
                restart = 0;
                // bottom-up
                for (i = 0; i < 5; i++) {
                    tempqij[i] = 0.0;
                    //printf("tempqij[%d] = %g\n", i, tempqij[i]);
                }
                int valL;
                int valL2;
                int comput[816] = {0};
                for (i = 0; i < 408; i++) {
                    for (j = 0; j < 6; j++) {
                        for (m = 0; m < 5; m++) {
                            if (m < j) {
                                valL = L[i][m]-1;
                                //tempqij[m] = qij[j][valL];
                                tempqij[m] = qij[comput[valL]][valL];
                                //printf("qij[%d][%d] = %g\n",comput[valL], valL, qij[comput[valL]][valL]);
                                //printf("%g \n", tempqij[m]);
                            } 
                            else if (m >= j) {
                                valL = L[i][m+1]-1;
                                //tempqij[m] = qij[j][valL];
                                tempqij[m] = qij[comput[valL]][valL];
                                //printf("qij[%d][%d] = %g\n",comput[valL], valL, qij[comput[valL]][valL]);
                                //printf("%g \n", tempqij[m]);
                            }
                        }
                        tempuij = tempqij[0];
                        double app;
                        double app1;
                        for(m = 1; m < 5; m++) {
                            //printf("tempuij = %g\n", tempuij);
                            app = sgn(tempuij) * sgn(tempqij[m]) * minabs(tempuij,tempqij[m]);
                            app1 = triangle(tempuij,tempqij[m]);
                            //printf("app = %g\n", app);
                            //printf("app1 = %g\n", app1);
                            //tempuij = ((sgn(tempuij) * sgn(tempqij[m]) * minabs(tempuij,tempqij[m])) + triangle(tempuij,tempqij[m]));
                            tempuij = app + app1;
                            //printf("tempqij[%d] = %g\n", m, tempqij[m]);
                            //printf("sgn(tempqij[%d]) = %d\n", m, sgn(tempqij[m]));
                            //printf("minabs(tempuij,tempqij[m]) = %g\n", minabs(tempuij,tempqij[m]));
                            //printf("triangle(tempuij,tempqij[m]) = %g\n", triangle(tempuij,tempqij[m]));
                            //sprintf("tempuij[%d] = %g\n", m, tempuij); 
                        }
                        //valL2 = L[i][j]-1;
                        //uij[i][valL2] = tempuij;
                        uij[i][j] = tempuij;
                        //printf("uij[%d][%d] = %g\n",i,j,uij[i][j]);
                    }
                    for (m = 0; m < 6; m++) {
                            comput[L[i][m] - 1] += 1;
                    }
                }

            // top-down
            for(i = 0; i < 3; i++) {
                temp1uij[i] = 0.0;
                //printf("temp1uij[%d] = %g\n", i, temp1uij[i]);
            }
            int comput1[408] = {0};
            for (j = 0; j < 816; j++) {
                for (i = 0; i < 3; i++) {
                    for (m = 0; m < 2; m++) {
                        if (m < i) { 
                            valL = M[j][m] - 1;
                            temp1uij[m] = uij[valL][comput1[valL]]; 
                        }
                        else if (m >= i) {
                            valL = M[j][m + 1] - 1;
                            temp1uij[m] = uij[valL][comput1[valL]];
                        }
                    }
                    temp1uij[2] = Lj[j];
                    //valL2 = M[j][i] - 1;
                    qij[i][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                }
                for (m = 0; m < 3; m++) {
                    comput1[M[j][m] - 1] += 1;
                }
            }

            // decision
            //printf("output = ");
            for (j = 0; j < 816; j++) {
                qj[j] = Lj[j] + uij[0][j] + uij[1][j] + uij[2][j]; 
                if (qj[j] >= 0) output[j] = 0;
                else if (qj[j] < 0) output[j] = 1;
                //printf("%d ", output[j]);
            }
            //printf("\n");

            // to check Hx=0
            //printf("checkbit = ");  
            //int compute;          
            for (i = 0; i < 408; i++) {
                checkbit[i] = 0;
                for (j = 0; j < 6; j++) {
                    //compute = H[i][j] * output[j];
                    //checkbit[i] += compute;
                    //checkbit[i] += (H[i][j] * output[j]);
                    checkbit[i] += output[L[i][j] - 1];
                }
                checkbit[i] = checkbit[i] % 2;
               // printf("%d",checkbit[i]);
            }
            //printf("\n");

            for (i = 0; i < 408; i++) {
                if (checkbit[i] == 0) restart += 1; // restart = 408 is success
            }
            stp = 0;
            if (k == 99 && restart != 408) {
                //printf("failure\n");
                stp = 1;
                s++;
            }
        }
        printf("k[%d] = %d\n", num, k);
        error = 0;
        for(i = 0; i < 816; i++) {
            if (output[i] != cod[i]) {
                error += 1;
                //printf("ouput[%d] = %d",i, output[i]);
            }
        }
        if (error != 0 && stp == 0) s++;
        restart = 0;
        //printf("error = %d\n", error);       
        totalerror += error;
        //if(num > 10000) break;
        //printf("totalerror[%d] = %d\n", num, totalerror);
    }
    double ber;
    ber = totalerror / (num * 816.0);
    printf("totalerror = %d\n", totalerror);
    printf("BER = %g\n", ber);
    bers[step] = ber;
    }


    for(step = 0; step < 6; step++) {
        printf("enb0s[%d] = %g\n",step, ebn0s[step]);
        printf("bers[%d] = %g\n",step, bers[step]);
    }

    FILE *outfp;
    
    outfp = fopen("resulttry.txt","w");
    for (i = 0; i < 6; i++) {
         fprintf(outfp,"%g ",ebn0s[i]);
         fprintf(outfp,"%g ",berscompare[i]);
         fprintf(outfp,"%g ",bers[i]);
         fprintf(outfp,"\n");
    }
    fclose(outfp);
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