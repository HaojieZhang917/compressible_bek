#include <stdio.h>
#include <stdlib.h>
#include<time.h>

int main()
{
    int number = 0;
    double sum = 0;
    int count = 0;
    printf("enter your number:");
    scanf("%d", &number);
    while (number != -1){
        sum += number;
        count++;
        printf("enter your number:");
        scanf("%d", &number);
    } 
    printf("the averange number of the number sequence is %f \n", sum/count);
    return 0;
}

