#include <stdio.h>
#include<stdlib.h>
#include <math.h>

int main() {
    double dB[6];
    double ber[6];
    double ber1[6];
    int i;

    dB[0] = 0.9844;     // Eb/N0 = 1.2544
    dB[1] = 1.289;      // Eb/N0 = 1.3456
    dB[2] = 1.584;      // Eb/N0 = 1.4401
    dB[3] = 1.868;
    dB[4] = 2.411;
    dB[5] = 3.046;

    ber[0] = 0.0583;
    ber[1] = 0.03052;
    ber[2] = 0.01011;
    ber[3] = 0.003;
    ber[4] = 8.34 * pow(10,-5);
    ber[5] = 2.799 * pow(10,-7);
    ber1[0] = 0.0456;
    ber1[1] = 0.03782;
    ber1[2] = 0.01111;
    ber1[3] = 0.00356;
    ber1[4] = 8.00 * pow(10,-5);
    ber1[5] = 2.5 * pow(10,-7);

    FILE *outfp;
    
    outfp = fopen("result.txt","w");
    for (i = 0; i < 6; i++) {
         fprintf(outfp,"%g ",dB[i]);
         fprintf(outfp,"%g ",ber[i]);
         fprintf(outfp,"%g ",ber1[i]);
         fprintf(outfp,"\n");
    }
    fclose(outfp);
    return 0;
}