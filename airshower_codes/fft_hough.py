import numpy as np
from numpy.fft import fft2, fftshift
import matplotlib.pyplot as plt
from PIL import Image, ImageOps
from skimage.transform import hough_line, hough_line_peaks

# Load and prepare the image
evtImage = Image.open('fardet_r00029307_s25_DDenergy_N23-10-24_v1_data.ultrashowerimages_ntuple_2090_yz.png')  # Replace with your path
evtImage = ImageOps.invert(evtImage.convert("L"))  # Invert grayscale

# Convert to array and crop center
array = np.asarray(evtImage, dtype=np.float32)
h, w = array.shape
halfH, halfW = h // 2, w // 2
r = min(halfH, halfW) - 5  # square crop radius
box = array[halfH - r:halfH + r, halfW - r:halfW + r]

# Compute FFT
fft_map = fftshift(fft2(box)) / (w * h)
fft_mag = np.log1p(np.square(1000 * np.abs(fft_map)))
#fft_mag = np.square(1000*np.abs(fft_map))
fft_mag[fft_mag < 8] = 0  # noise threshold

# Hough transform
hspace, angles, dists = hough_line(fft_mag, theta=np.linspace(-np.pi/2, np.pi/2, 3600))
peaks, peakAngles, peakDists = hough_line_peaks(hspace, angles, dists, num_peaks=10, min_angle=1, min_distance=1)

# --- Plotting ---
fig, axs = plt.subplots(1, 2, figsize=(12, 7))

# 1. Original (inverted) image box
#axs[0].imshow(box, cmap='gray')
#axs[0].set_title('Cropped Center (Inverted)')
#axs[0].axis('off')

# 2. FFT Map
axs[0].imshow(fft_mag, cmap='inferno')
#axs[1].set_title('FFT Log Magnitude')
axs[0].axis('off')

# 3. Hough lines on FFT Map
axs[1].imshow(fft_mag, cmap='gray')
axs[1].axis('off')

# Lock axis limits to prevent resizing from plotted lines
axs[1].set_xlim(0, fft_mag.shape[1])
axs[1].set_ylim(fft_mag.shape[0], 0)  # flip Y-axis for image coordinates

# Overlay Hough lines
for angle, dist in zip(peakAngles, peakDists):
    x = np.linspace(0, fft_mag.shape[1], 1000)
    y = (dist - x * np.cos(angle)) / np.sin(angle)
    axs[1].plot(x, y, '-r', linewidth=1)

plt.tight_layout()
plt.savefig("fft_hough_yz.png")
plt.savefig("fft_hough_yz.eps")
plt.savefig("fft_hough_yz.pdf")
#plt.show()
