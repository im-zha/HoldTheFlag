# Vòng 1 — Mệnh lệnh xuất kích
## Given
File challenge cung cấp một bài toán **Discrete Logarithm Problem (DLP)** trên nhóm số nguyên modulo `p`:
```
The flag has been encoded as a secret exponent x, where:

  h = g^x mod p

Your job: find x. Convert it from integer to bytes to get the flag.

p = 2247297901375864750461918215908397462554622316861885697469812609208482671796519989041708819088327406276727713325926984783654453580598954278480933663642628407743971547337698858815665335450301828703735517279088545252493204317839245564990751219121256870972082020870742841156981131437941510099651935803198999627
g = 797032223149531285607971355158192150926555680763341545661871333078567437408820428834941794887222745425244137303956776516613894961185972456200978813725184212858037439841778740577859187417009007458444768864113847557627519817788214445986970343968064264731545268364867143251499664214088696995903599781563984116
h = 468236306785559227242545208642961660487424736475509451324848496240713063884248400928453707063093512040762457493692466389598131256501402625717666645351148316070083773749717409678738689693550470750839264535543953759648092479282888291678698194889269011735062078161680411521433536665807130299043504305373037050
```

## Goal
Tìm số nguyên bí mật `x` thỏa mãn phương trình trên, sau đó chuyển đổi `x` từ số nguyên sang bytes để thu được flag theo format `FLAG{[A-Z0-9_]+}`.

## Solution
- **Bước 1 — Nhận diện lỗ hổng: Smooth `p-1`**

    Ta thấy `p - 1` phân tích thành 45 thừa số nguyên tố, tất cả đều nhỏ hơn 2²⁴ ≈ 16.7 triệu.

    Đây gọi là **B-smooth number**.

    => Sai lầm: trong hệ thống DLP an toàn (như Diffie-Hellman chuẩn), `p` phải là safe prime (tức $p = 2q + 1$ với $q$ nguyên tố lớn), để `p - 1` có ít nhất một thừa số nguyên tố rất lớn. Khi `p - 1` smooth, thuật toán Pohlig-Hellman có thể khai thác.

- **Bước 2 — Pohlig-Hellman Attack**

    **Ý tưởng cốt lõi của Pohlig-Hellman:** thay vì giải DLP trên nhóm có order `p - 1` khổng lồ, chia bài toán thành các DLP nhỏ trên các nhóm con có order là từng thừa số nguyên tố $q_i$, sau đó ghép lại bằng Định lý số dư Trung Hoa (CRT):

    $$x \equiv x_i \pmod{q_i} \quad \text{với mỗi thừa số } q_i \mid (p - 1)$$

    Vì tất cả $q_i < 2^{24}$, mỗi DLP con có thể giải bằng \textbf{Baby-step Giant-step} trong $O(\sqrt{q_i}) \approx O(2^{12})$ bước — cực kỳ nhanh.

- **Bước 3 - Decode flag**

    Sau khi Pohlig-Hellman tìm ra `x` (là một số nguyên lớn), ta cần chuyển số đó thành chuỗi ký tự.

## Code
```py
import math
from functools import reduce


p = 2247297901375864750461918215908397462554622316861885697469812609208482671796519989041708819088327406276727713325926984783654453580598954278480933663642628407743971547337698858815665335450301828703735517279088545252493204317839245564990751219121256870972082020870742841156981131437941510099651935803198999627
g = 797032223149531285607971355158192150926555680763341545661871333078567437408820428834941794887222745425244137303956776516613894961185972456200978813725184212858037439841778740577859187417009007458444768864113847557627519817788214445986970343968064264731545268364867143251499664214088696995903599781563984116
h = 468236306785559227242545208642961660487424736475509451324848496240713063884248400928453707063093512040762457493692466389598131256501402625717666645351148316070083773749717409678738689693550470750839264535543953759648092479282888291678698194889269011735062078161680411521433536665807130299043504305373037050


def bsgs(g, h, p, order):
    m = math.isqrt(order) + 1
    table = {}
    gj = 1
    for j in range(m):
        table[gj] = j
        gj = gj * g % p
    g_inv_m = pow(g, p - 1 - m, p)
    gamma = h
    for i in range(m):
        if gamma in table:
            x = i * m + table[gamma]
            if x < order:
                return x % order
        gamma = gamma * g_inv_m % p
    return None


def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1 = extended_gcd(b % a, a)
    return gcd, y1 - (b // a) * x1, x1


def crt(residues, moduli):
    M = reduce(lambda a, b: a * b, moduli)
    x = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        _, inv, _ = extended_gcd(Mi % m, m)
        x += r * Mi * inv
    return x % M


def pohlig_hellman(g, h, p, factors):
    n = p - 1
    residues, moduli = [], []
    for q, e in factors.items():
        qe = q ** e
        gi = pow(g, n // qe, p)
        hi = pow(h, n // qe, p)
        if e == 1:
            xi = bsgs(gi, hi, p, q) or 0
        else:
            xi = 0
            gamma = pow(gi, qe // q, p)
            for k in range(e):
                hk = pow(pow(g, xi, p) * pow(h, p - 2, p) % p, n // q**(k+1), p)
                dk = bsgs(gamma, pow(hk, p - 2, p), p, q) or 0
                xi = (xi + dk * q**k) % qe
        residues.append(xi)
        moduli.append(qe)
    return crt(residues, moduli)


# Bước 1: Phân tích p-1
from sympy import factorint
factors = factorint(p - 1)
print(len(factors))

# Bước 2: Giải DLP
x = pohlig_hellman(g, h, p, factors)

# Bước 3: Verify & decode
assert pow(g, x, p) == h, "Verification failed!"
flag = x.to_bytes((x.bit_length() + 7) // 8, "big").decode("utf-8")
print(flag)
```

## Flag
`FLAG{HQ604}`

# VÒNG 2 — Tọa độ điểm mù
## Given
File ảnh `holdtheflag.png` chỉ là một màu đen, đi kèm với đoạn mô tả chứa manh mối quan trọng: "Ở đây tối quá, chẳng nhìn thấy gì cả."

## Goal
Phân tích bức ảnh, tìm cách làm lộ ra các thông tin bị giấu trong vùng tối để tìm được tọa độ/tên địa danh, từ đó định dạng lại thành Flag hợp lệ theo chuẩn `FLAG{[A-Z0-9_]+}`.

## Solution
- Ảnh tối nhưng các text label vẫn đọc được khi tăng độ sáng/contrast. 

- Ta đọc 3 tên theo thứ tự từ nổi bật nhất -> ít nổi bật hơn (theo sự kiện lịch sử: Gạc Ma là nơi xảy ra trận chiến chính):

    ```
    GAC_MA -> CO_LIN -> LEN_DAO
    ```

## Flag
`FLAG{GAC_MA_CO_LIN_LEN_DAO}`