#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *input_file1 = fopen("test_data.txt", "r");
    int num;
    int input_i = 0;
    while(fscanf(input_file1, "%d", &num) == 1) {
        // input_one[input_i] = num;
        if (input_i == 0) {
            printf("%d\n", num);
        }
        input_i += 1;
    }




}

