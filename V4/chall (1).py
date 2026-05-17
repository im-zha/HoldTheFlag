import os
from Crypto.Cipher import AES

KEY = ?
class AES_CTR:
    def __init__(self, key, step_up=False):
        self.key = key
        self.cipher_core = AES.new(self.key, AES.MODE_ECB)
        
        self.value = os.urandom(16).hex()
        self.step = 1
        self.stup = step_up

    def increment(self):
        if self.stup:
            self.newIV = hex(int(self.value, 16) + self.step)
        else:
            self.newIV = hex(int(self.value, 16) - self.stup) 
            
        self.value = self.newIV[2:len(self.newIV)]
        counter_block = bytes.fromhex(self.value.zfill(32))
        
        return self.cipher_core.encrypt(counter_block)

    def encrypt(self, data: bytes) -> bytes:
        out = bytearray()
        
        for i in range(0, len(data), 16):
            block = data[i:i+16]
            keystream = self.increment()
            
            xored = bytes(a ^ b for a, b in zip(block, keystream))
            out.extend(xored)
            
        return bytes(out)


def encrypt_challenge():

    cipher = AES_CTR(KEY, step_up=False)

    with open("flag.png", 'rb') as f:
        plaintext = f.read()

    ciphertext = cipher.encrypt(plaintext)

    encrypted_hex = ciphertext.hex()  
    with open("output.txt", 'w') as f:
        f.write(encrypted_hex)

encrypt_challenge()