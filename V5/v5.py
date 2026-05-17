# Chuỗi hex tìm được trong header của file
hex_string = "0b19636f651d08737e6b1c197f757a0f0e7a79710616"
ciphertext = bytes.fromhex(hex_string)

# Khóa giải mã chính là tên file
key = b"HQ604"

# Tiến hành giải mã XOR
plaintext = bytearray()
for i in range(len(ciphertext)):
    # XOR từng byte của ciphertext với byte tương ứng của key (lặp vòng)
    decrypted_byte = ciphertext[i] ^ key[i % len(key)]
    plaintext.append(decrypted_byte)

# In kết quả
flag_content = plaintext.decode('utf-8')
print(f"FLAG{{{flag_content}}}")