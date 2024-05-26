#include <string>
#include <iostream>
#include <algorithm>
#include "ccipher.h"


// -------------------------------------------------------

string rot_alpha(int offset) {
    string alpha = "abcdefghijklmnopqrstuvwxyz";
    rotate_string(alpha, offset);
    return alpha;
}

// create new sub class and pass it the rotated string
// calling Cipher outside body bc parent is built first and needs to
// be given it's argument before child class body is run
CCipher::CCipher(int offset) : Cipher(rot_alpha(offset)) {}


// ./cipher -m c -i test2/caesar/caesar.txt -o test2/caesar/caesar_liam_out.txt
// -------------------------------------------------------


// Rotates the input string in_str by rot positions
void rotate_string(string& in_str, int rot)
{
    // TODO: You will likely need this function. Implement it.
    string firstHalf = "";
    string secondHalf = "";
    int normalizedOffset = rot % 26;
    int counter = 0;
    for (char c : in_str) {
        if (counter < normalizedOffset) {
            firstHalf += c;
        } else {
            secondHalf += c;
        }
        counter += 1;
    }
    reverse(firstHalf.begin(), firstHalf.end());
    reverse(secondHalf.begin(), secondHalf.end());
    string tmp = firstHalf + secondHalf;
    reverse(tmp.begin(), tmp.end());
    in_str = tmp;
}
