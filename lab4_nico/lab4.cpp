#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double epsilon = 1e-10;

void thomasAlgorithm(double *a, double *b, double *c, double *d, double *x, int n)
{
    if(n == 1){
        x[0] = d[0]/b[0];
    } else {
        c[0] = c[0]/b[0];
        d[0] = d[0]/b[0];
        for(int i = 1; i < n; i++){
            c[i] = c[i]/(b[i] - a[i]*c[i-1]);
            d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
        }
        x[n-1] = d[n-1];
        for(int i = n-2; i >= 0; i--)
            x[i] = d[i] - c[i]*x[i+1];
    }
}

bool drum1(double x, double y) // return true for inner region
{
    if((y > 0)&&(y < 1/3.0)&&(x > 2/3.0 + epsilon)&&(y > x - 2/3.0 + epsilon)) return true;
    if((y >= 1/3.0)&&(y < 2/3.0)&&(y > -x + 2/3.0 + epsilon)&&(x < 1 - epsilon)) return true;
    if((y >= 2/3.0)&&(y < 1)&&(y < x + 2/3.0 - epsilon)&&(x < 1/3.0 - epsilon)) return true;
    return false;
}
bool drum2(double x, double y) // return true for inner region
{

    if((y > 0)&&(y < 1/3.0)&&(y > -x + 2/3.0 + epsilon)&&(x < 2/3.0 - epsilon)) return true;
    if((y >= 1/3.0)&&(y < 2/3.0)&&(y > -x + 2/3.0 + epsilon)&&(y < -x + 4/3.0 - epsilon)) return true;
    if((y >= 2/3.0)&&(y < 1)&&(x > 0)&&(x < 1/3.0 - epsilon)) return true;
    return false;
}

bool drum3(double x, double y) // return true for inner region
{
    if(((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) < 0.24)&&((y < 1/2.0 - epsilon)||(x < 1/3.0 - epsilon))) return true;
    return false;
}