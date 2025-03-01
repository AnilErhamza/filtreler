#include <math.h>
#include <stdio.h>


// kalmanın kullandığı fonksiyonlar ve kaydedilmesi gereken değişkenler burda tanımlı

struct array3{
	float a0[3];
	float a1[3];
	float a2[3];
}array3;

struct array{
	float b0;
	float b1;
	float b2;
}array;

struct array3 H = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

struct array3 Q = {{0.1, 0, 0}, {0, 0.1, 0}, {0, 0, 0.1}}; // process noise covariance matrix

struct array3 R = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}}; // measurement noise covariance matrix

struct array xx; 

struct array3 P = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

struct array3 A = {{1,0,0},{0,1,0},{0,0,1}};



struct array3 multiplyMatrix(struct array3 m1, struct array3 m2)
{
	struct array3 result;
	result.a0[0]=m1.a0[0]*m2.a0[0] + m1.a0[1]*m2.a1[0] + m1.a0[2]*m2.a2[0];
	result.a0[1]=m1.a0[0]*m2.a0[1] + m1.a0[1]*m2.a1[1] + m1.a0[2]*m2.a2[1];
	result.a0[2]=m1.a0[0]*m2.a0[2] + m1.a0[1]*m2.a1[2] + m1.a0[2]*m2.a2[2];
	
	result.a1[0]=m1.a1[0]*m2.a0[0] + m1.a1[1]*m2.a1[0] + m1.a1[2]*m2.a2[0];
	result.a1[1]=m1.a1[0]*m2.a0[1] + m1.a1[1]*m2.a1[1] + m1.a1[2]*m2.a2[1];
	result.a1[2]=m1.a1[0]*m2.a0[2] + m1.a1[1]*m2.a1[2] + m1.a1[2]*m2.a2[2];
	
	result.a2[0]=m1.a2[0]*m2.a0[0] + m1.a2[1]*m2.a1[0] + m1.a2[2]*m2.a2[0];
	result.a2[1]=m1.a2[0]*m2.a0[1] + m1.a2[1]*m2.a1[1] + m1.a2[2]*m2.a2[1];
	result.a2[2]=m1.a2[0]*m2.a0[2] + m1.a2[1]*m2.a1[2] + m1.a2[2]*m2.a2[2];
	
	return result;	
}

struct array multiplyMatrix2(struct array3 m1, struct array m2)
{
	struct array result;
	result.b0= m1.a0[0]*m2.b0 + m1.a0[1]*m2.b0 + m1.a0[2]*m2.b0;
	result.b1= m1.a1[0]*m2.b1 + m1.a1[1]*m2.b1 + m1.a1[2]*m2.b1;
	result.b2= m1.a2[0]*m2.b2 + m1.a2[1]*m2.b2 + m1.a2[2]*m2.b2;
	return result;
}

struct array3 transpose2(struct array3 A)
{
	struct array3 B;
   B.a0[0] = A.a0[0];
   B.a0[1] = A.a1[0];
   B.a0[2] = A.a2[0];
   
   B.a1[0] = A.a0[1];
   B.a1[1] = A.a1[1];
   B.a1[2] = A.a2[1];
   
   B.a2[0] = A.a0[2];
   B.a2[1] = A.a1[2];
   B.a2[2] = A.a2[2];
   
   return B;
}

struct array3 sum(struct array3 a, struct array3 b)
{
  struct array3 sum;
  
  	sum.a0[0] = a.a0[0] + b.a0[0];
  	sum.a0[1] = a.a0[1] + b.a0[1];
 	sum.a0[2] = a.a0[2] + b.a0[2];
 	
	sum.a1[0] = a.a1[0] + b.a1[0];
  	sum.a1[1] = a.a1[1] + b.a1[1];
  	sum.a1[2] = a.a1[2] + b.a1[2];

	sum.a2[0] = a.a2[0] + b.a2[0];
  	sum.a2[1] = a.a2[1] + b.a2[1];
  	sum.a2[2] = a.a2[2] + b.a2[2];
  
  	return sum;
}

struct array sum2(struct array a, struct array b)
{
  struct array sum;
  
  	sum.b0 = a.b0 + b.b0;
  	sum.b1 = a.b1 + b.b1;
 	sum.b2 = a.b2 + b.b2;
 	
  	return sum;
}

struct array3 sub(struct array3 a, struct array3 b)
{
  struct array3 sub;
  
  	sub.a0[0] = a.a0[0] - b.a0[0];
  	sub.a0[1] = a.a0[1] - b.a0[1];
 	sub.a0[2] = a.a0[2] - b.a0[2];
 	
	sub.a1[0] = a.a1[0] - b.a1[0];
  	sub.a1[1] = a.a1[1] - b.a1[1];
  	sub.a1[2] = a.a1[2] - b.a1[2];

	sub.a2[0] = a.a2[0] - b.a2[0];
  	sub.a2[1] = a.a2[1] - b.a2[1];
  	sub.a2[2] = a.a2[2] - b.a2[2];
  
  	return sub;
}

struct array sub2(struct array a, struct array b)
{
  struct array sub;
  
  	sub.b0 = a.b0 - b.b0;
  	sub.b1 = a.b1 - b.b1;
 	sub.b2 = a.b2 - b.b2;
 	
  	return sub;
}

float determinant(float a[3][3], float k)
{
  float s = 1, det = 0, b[3][3];
  int i, j, m, n, c;
  if (k == 1)
  {
    return (a[0][0]);
  }
  else
  {
    det = 0;
    for (c = 0; c < k; c++)
    {
      m = 0;
      n = 0;
      for (i = 0; i < k; i++)
      {
        for (j = 0; j < k; j++)
        {
          b[i][j] = 0;
          if (i != 0 && j != c)
          {
            b[m][n] = a[i][j];
            if (n < (k - 2))
              n++;
            else
            {
              n = 0;
              m++;
            }
          }
        }
      }
      det = det + s * (a[0][c] * determinant(b, k - 1));
      s = -1 * s;
    }
  }

  return det;
}

struct array3 transpose(float num[3][3], float fac[3][3], float r)
{
	float inverse[3][3];
  int i, j;
  float b[3][3], d;

  for (i = 0; i < r; i++)
  {
    for (j = 0; j < r; j++)
    {
      b[i][j] = fac[j][i];
    }
  }
  d = determinant(num, r);
  for (i = 0; i < r; i++)
  {
    for (j = 0; j < r; j++)
    {
      inverse[i][j] = b[i][j] / d;
    }
  }
  
  	struct array3 inverse2;
  
  	inverse2.a0[0] = inverse[0][0]   ;
	inverse2.a0[1] = inverse[1][0]   ;
	inverse2.a0[2] = inverse[2][0]  ;
	
	inverse2.a1[0] = inverse[0][1]   ;
	inverse2.a1[1] = inverse[1][1]   ;
	inverse2.a1[2] = inverse[2][1]   ;
	
	inverse2.a2[0] = inverse[0][2]   ;
	inverse2.a2[1] = inverse[1][2]  ;
	inverse2.a2[2] = inverse[2][2]  ;
  
  	return inverse2;
  
}

struct array3 cofactor(struct array3 num)
{
	struct array3 inv;
	float num2[3][3];
	
	num2[0][0] = num.a0[0] ;
	num2[0][1] = num.a1[0] ;
	num2[0][2] = num.a2[0] ;
	
	num2[1][0] = num.a0[1] ;
	num2[1][1] = num.a1[1] ;
	num2[1][2] = num.a2[1] ;
	
	num2[2][0] = num.a0[2] ;
	num2[2][1] = num.a1[2] ;
	num2[2][2] = num.a2[2] ;
	
	
  float b[3][3], fac[3][3];
  int p, q, m, n, i, j, f = 3;
  for (q = 0; q < f; q++)
  {
    for (p = 0; p < f; p++)
    {
      m = 0;
      n = 0;
      for (i = 0; i < f; i++)
      {
        for (j = 0; j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num2[i][j];
            if (n < (f - 2))
              n++;
            else
            {
              n = 0;
              m++;
            }
          }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  inv = transpose(num2, fac, 3);
  
  return inv;
}

struct array3 makenull(struct array3 decoy1){
	decoy1.a0[0] = 0;
	decoy1.a0[1] = 0;
	decoy1.a0[2] = 0;
	
	decoy1.a1[0] = 0;
	decoy1.a1[1] = 0;
	decoy1.a1[2] = 0;
	
	decoy1.a2[0] = 0;
	decoy1.a2[1] = 0;
	decoy1.a2[2] = 0;
	
	return decoy1;
}

/*Finding transpose of matrix*/
