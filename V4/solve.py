# PNG Header (16 bytes đầu tiên của mọi file PNG)
png_header = bytes.fromhex("89504E470D0A1A0A0000000D49484452")

# Đọc dữ liệu đã mã hóa
with open(r"D:\HoldTheFlag\V4\output.txt", "r") as f:
    ciphertext_hex = f.read().strip()
ciphertext = bytes.fromhex(ciphertext_hex)

# Lấy 16 byte ciphertext đầu tiên
first_block = ciphertext[:16]

# Tìm Keystream (vì Keystream bị cố định)
# Keystream = Ciphertext ^ Plaintext
keystream = bytes(a ^ b for a, b in zip(first_block, png_header))

# Khôi phục file ảnh (Giải mã)
decrypted_data = bytearray()
for i in range(0, len(ciphertext), 16):
    block = ciphertext[i:i+16]
    current_keystream = keystream[:len(block)]
    decrypted_block = bytes(a ^ b for a, b in zip(block, current_keystream))
    decrypted_data.extend(decrypted_block)

# Lưu kết quả
with open(r"D:\HoldTheFlag\V4\flag_decrypted.png", "wb") as f:
    f.write(decrypted_data)