#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include "kcipher.h"

/* Helper function definitions
 */

int lengthStr(string word) {
    int count = 0;
    for (char c : word) {
        if (c != ' ') {
            count += 1;
        }
    }
    return count;
}

struct KCipher::KCipherCheshire {
    vector<string> run_keys;
    int page;
    int page_count = 0;
};

void KCipher::set_id(int id) {
    ksmile->page = id;
}

void KCipher::add_key(string word) {
    if (lengthStr(word) == 0){
        exit(EXIT_FAILURE);
    }
    (ksmile->run_keys).push_back(word);
    ksmile->page_count += 1;
}

KCipher::KCipher(string word) {
    if (lengthStr(word) == 0){
        exit(EXIT_FAILURE);
    }
    ksmile = new KCipherCheshire;
    (ksmile->run_keys).push_back(word);
    ksmile->page_count += 1;
}

KCipher::KCipher() {
    ksmile = new KCipherCheshire;
    string allAs(MAX_LENGTH, 'a');
    ksmile->page = 0;
    (ksmile->run_keys).push_back(allAs);
}


// encrypt implementation
// Done!
string KCipher::encrypt(string plain_text) {
    int key_i = 0;
    string base_alpha = "abcdefghijklmnopqrstuvwxyz";
    string rot_alpha = "abcdefghijklmnopqrstuvwxyz";
    string running_key = (ksmile->run_keys)[ksmile->page];
    string ret = "";
    if (lengthStr(running_key) < lengthStr(plain_text)) {
        exit(EXIT_FAILURE);
    }

    for (char c : plain_text) {
        if (c == ' ') {
            ret += " ";
        }
        else {
            int row = find_pos(base_alpha, c);
            rotate_string(rot_alpha, row);
            if (running_key[key_i] == ' ') {
                key_i += 1;
            }
            int col = find_pos(base_alpha, running_key[key_i]);
            if (isupper(c)) {
                ret += toupper(rot_alpha[col]);
            } else{
                ret += rot_alpha[col];
            }
            rot_alpha = base_alpha;
            key_i += 1;
        }
    }
    return ret;
}
// ./cipher -m k -i test2/running_key/running_key.txt -o test2/running_key/running_key_liam_out.txt
// set args -m k -i test2/running_key/running_key.txt -o test2/running_key/running_key_liam_out.txt



KCipher::~KCipher() {
    delete ksmile;
}

// decrypt implementation
string KCipher::decrypt(string cipher_text) {
    int key_i = 0;
    string base_alpha = "abcdefghijklmnopqrstuvwxyz";
    string rot_alpha = "abcdefghijklmnopqrstuvwxyz";
    string running_key = (ksmile->run_keys)[ksmile->page];
    string retStr = "";
    for (char c : cipher_text) {
        if (c == ' ') {
            retStr += " ";
        } else {
            // where k should be in the alphabet
            if (running_key[key_i] == ' ') {
                key_i += 1;
            }
            char r_letter = running_key[key_i];
            int offset = r_letter - tolower(c);
            // where k is in the alphabet
            rotate_string(rot_alpha, (26 - (offset)%26));
            if (isupper(c)) {
                retStr += toupper(rot_alpha[0]);
            } else {
                retStr += rot_alpha[0];
            }
            rot_alpha = base_alpha;
            key_i += 1;
        }
    }
    return retStr;
}
// -------------------------------------------------------
// Running Key Cipher implementation

// -------------------------------------------------------

