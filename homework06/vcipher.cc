#include <string>
#include <iostream>
#include <vector>
#include "kcipher.h"
#include "vcipher.h"


// -------------------------------------------------------
// Running Key Cipher implementation
// -------------------------------------------------------
string multStr(string word) {
    string added = "";
    for (int i = 0; i < 30; i++) {
        added += word;
    }
    return added;
}

VCipher::VCipher(string word) : KCipher(multStr(word)) {}
VCipher::VCipher() : KCipher() {}

// ./cipher -m v -i test2/vigenere/vigenere.txt -o test2/vigenere/vigenere_liam_out.txt
// set args -m v -i test2/vigenere/vigenere.txt -o test2/vigenere/vigenere_liam_out.txt
