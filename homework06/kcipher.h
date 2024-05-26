#ifndef KCIPHER_H_
#define KCIPHER_H_
#include "cipher.h"
#include "ccipher.h"

using namespace std;

const unsigned int MAX_LENGTH = 100;

/* A class that implements a
   Running key cipher class. It 
   inherts class Cipher */
// TODO: Implement this function
class KCipher : public Cipher {
    private:
        struct KCipherCheshire;
        KCipherCheshire *ksmile;
    public:
        KCipher();
        KCipher(string word);
        ~KCipher();
        void add_key(string word);
        void set_id(int id);
        string encrypt(string plain_text) override;
        string decrypt(string cipher_text) override;
};

#endif