#include<stdio.h>

double a,b;
int main(){
	while(scanf("%lf %lf",&a,&b) > 0){
		printf("%lf\n",b-a);
	}
}
