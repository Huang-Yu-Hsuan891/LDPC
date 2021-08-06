#include <stdio.h>
#include <math.h>

int sgn(double L);
double minabs(double a, double b);
double triangle(double L1, double L2);

int main() {
    //int x = sgn(-1);
    //double x = triangle(0.2, -0.1);
    double x;
    x = sgn(2.826) * sgn(-0.79049) * minabs(2.826,-0.79049) + triangle(2.826,-0.79049);
    printf("x = %g\n",x);
    return 0;
}


int sgn (double L){
    if (L >= 0) return 1;
    else return -1;
}

double minabs(double a, double b) {
    if(a <= 0) a = (-1) * a;
    else a = a;
    if(b <= 0) b = (-1) * b;
    else b = b;
    if(a>=b) return b;
    else return a;
}

double triangle(double L1, double L2) {
    double temp1, temp2;
    double ope1,ope2;
    double answer;
    ope1 = L1 + L2;
    ope2 = L1 - L2;
    printf("ope1 = %g; ope2 = %g\n", ope1,ope2);
    if (ope1 <= 0) ope1 = ope1;
    else ope1 = (-1) * ope1;
    if (ope2 <= 0) ope2 = ope2;
    else ope2 = (-1) * ope2;
    printf("ope1 = %g; ope2 = %g\n", ope1,ope2);    
    temp1 = 1 + exp(ope1);
    temp2 = 1 + exp(ope2);
    printf("temp1 = %g\n", temp1);
    printf("temp2 = %g\n", temp2);
    printf("temp1/temp2 = %g\n", temp1/temp2);

    answer = log(temp1 / temp2);
    printf("answer = %g\n", answer);
    return answer;
}