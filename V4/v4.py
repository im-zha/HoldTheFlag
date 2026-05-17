import matplotlib.pyplot as plt
from scipy.io import wavfile

# Đọc file audio
sample_rate, data = wavfile.read(r"D:\HoldTheFlag\V4\challenge.wav")

# Nếu audio có 2 kênh (stereo), chỉ lấy 1 kênh để xử lý
if len(data.shape) > 1:
    data = data[:, 0]

# Vẽ và cấu hình Spectrogram
plt.figure(figsize=(12, 6))
# NFFT càng cao thì độ phân giải dải tần càng chi tiết
plt.specgram(data, Fs=sample_rate, NFFT=1024, cmap='magma') 

plt.title('Hidden Message in Spectrogram')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.ylim(0, 15000) # Giới hạn tần số để text hiện ra rõ nhất (có thể tinh chỉnh)

# Lưu kết quả ra file ảnh để đọc Flag
plt.savefig(r"D:\HoldTheFlag\V4\spectrogram_flag.png")
print("Đã giải mã âm thanh! Hãy mở file spectrogram_flag.png để lấy cờ.")