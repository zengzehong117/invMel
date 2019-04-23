#include<stdio.h>
#include<complex.h>
#include<math.h>
#include<stdlib.h>

void scan(double complex *a, int n)
{
    for(int i = 0; i < n;i++)
    {
        double real,imag;
        scanf("%lf %lf",&real,&imag);
        a[i] = real + imag*I;
    }
}
void print(double complex *a, int n,int key)
{
    for(int i = 0; i < n;i++)
    {
        if(key == 1)
            a[i] = a[i]/n;
        printf("%.6lf %.6lf\n",creal(a[i]),cimag(a[i]));
    }

}
double complex * FFT(double complex *a, int length,int key)
{
    int l = length;
    if(l == 1)
        return a;
    double pi = acos(-1);
    double complex w_l= cos(2*pi/l) + sin(2*pi/l)*I;
    if(key == 1)
        w_l = 1/w_l;
    double complex w = 1.0;
    double complex *a_even = (double complex *) malloc((sizeof(double complex))*l/2);
    double complex *a_odd = (double complex *) malloc((sizeof(double complex))*l/2);
    for(int k = 0; k<= l/2 - 1;++k)
    {
        a_even[k] = a[2*k];
        a_odd[k] = a[2*k+1];
    }
    double complex *y_even = FFT(a_even,l/2,key);
    double complex *y_odd = FFT(a_odd,l/2,key);
    double complex *y = (double complex *) malloc((sizeof(double complex))*l);
    for(int k = 0; k <= l/2-1;++k)
    {
        y[k] = y_even[k] + w*(y_odd[k]);
        y[k+l/2] = y_even[k] - w*(y_odd[k]);
        w = w*w_l;
    }
    free(a_even);
    free(a_odd);

    return y;
}



/*
int main()
{
    int t;
    scanf("%d",&t);
    for(int i = 0; i<t;++i)
    {
        int key,n;
        scanf("%d %d",&key,&n);
        double complex *y = (double complex *) malloc(sizeof(double complex)*n);
        double complex *a = (double complex *) malloc(sizeof(double complex)*n);

        scan(a,n);
        y = FFT(a,n,key);
        print(y,n,key);

    }
    return 0;
}
*/
