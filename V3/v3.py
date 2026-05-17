import base64

encoded_str = "SU9ESntXVURRX1lEUV9TS1hSUUp9"
decoded_b64 = base64.b64decode(encoded_str).decode('utf-8')

flag = ""
for char in decoded_b64:
    if char.isalpha():
        shift_base = 65 if char.isupper() else 97
        flag += chr((ord(char) - shift_base - 3) % 26 + shift_base)
    else:
        flag += char

print(flag)