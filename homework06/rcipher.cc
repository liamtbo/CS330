#include <string>
#include <iostream>
#include "rcipher.h"

// -------------------------------------------------------
// ROT13 Cipher implementation
// -------------------------------------------------------

// create new sub class and pass it the rotated string
RCipher::RCipher() : CCipher(13) {}

// ./cipher -m r -i test2/rot13/rot13.txt -o test2/rot13/rot13_liam_out.txt
