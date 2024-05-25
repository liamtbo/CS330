#include <stdio.h>
#include <string.h>

int main() {
    char *tmp1 = "butter";
    char *tmp2 = "butter";
    char *tmp3 = strcat(tmp1, tmp2);
    printf("%s\n", tmp1);
    printf("%s\n", tmp2);
    printf("%s\n", tmp3);

}