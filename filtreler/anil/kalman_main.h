
#include "kalman.h"
#include <stdio.h>

// kalman algoritması, 3 float alıp array veriyo

// burda array dediğim kalman.h da oluşturduğum struct

//aşşağıda printler var silebilirsin


struct array kalman(float yaw, float pitch, float roll)
{
    struct array z = {yaw, pitch, roll};
    struct array xp;
    struct array3 Pp;
    struct array3 At;
    struct array3 Ht;
    struct array3 K;
    struct array3 decoy1, decoy2, decoy3, decoy4;
    struct array decoy5, decoy6, decoy7, decoy8;
  

	xp = multiplyMatrix2(A, xx);
	  
    decoy1 = multiplyMatrix(A, P);

    At = transpose2(A);

    decoy2 = multiplyMatrix(decoy1, At);

   	Pp = sum(decoy2, Q);
	

    
    //
    
    decoy1 = makenull(decoy1);
	decoy2 = makenull(decoy2);
    
	// 
	
    decoy1 = multiplyMatrix(H, Pp);

    Ht = transpose2(H);

    decoy2 = multiplyMatrix(decoy1, Ht);

    decoy3 = sum(decoy2, R);

    decoy4 = cofactor(decoy3);

    decoy1 = multiplyMatrix(Pp, Ht);

    K = multiplyMatrix(decoy1, decoy4);

    //
    
    decoy1 = makenull(decoy1);
    decoy2 = makenull(decoy2);
    decoy3 = makenull(decoy3);
    decoy4 = makenull(decoy4);
    
    //

    decoy5 = multiplyMatrix2(H, xp);

    decoy6 = sub2(z, decoy5);

    decoy7 = multiplyMatrix2(K, decoy6);

    xx = sum2(xp, decoy7);

    //
    
    decoy1 = makenull(decoy1);
    decoy2 = makenull(decoy2);
    decoy3 = makenull(decoy3);
    decoy4 = makenull(decoy4);
    
    //

    decoy1 = multiplyMatrix(K, H);

    decoy2 = multiplyMatrix(decoy1, Pp);

    P = sub(Pp, decoy2);

    // printler x i yazdırıyo değiştikten sonra test için sadece
    printf("x=%f\n", xx.b0);
    printf("x=%f\n", xx.b1);
    printf("x=%f\n", xx.b2);
    printf("\n");
    
    return xx;
}
