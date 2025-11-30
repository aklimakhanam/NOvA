import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# define data points
data = {
    'PID': [0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6],
    'S_over_sqrtB_et_fc': [3.953, 4.014, 4.096, 4.213, 4.328, 4.459, 4.451, 4.487, 4.772, 4.646, 4.330],
    'S_over_sqrtB_fc': [3.423, 3.494, 3.595, 3.696, 3.771, 3.992, 4.013, 4.050, 4.223, 4.090, 3.80]
}

# Create a pandas DataFrame
df = pd.DataFrame(data)

min_x = df['PID'].min()
max_x = df['PID'].max()

# plot
plt.figure(figsize=(8, 6))
#plt.plot(df['PID'], df['S_over_sqrtB_et_fc'], marker='o', linestyle='-', color='blue', label='ET FC')
plt.plot(df['PID'], df['S_over_sqrtB_fc'], marker='s', linestyle='-', color='red')
plt.title(r'FOM vs PID')
plt.xlabel('PID')
plt.ylabel(r'$S / \sqrt{B}$')
plt.xticks(np.arange(min_x, max_x + 0.02, 0.02))
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("fom_vs_pid.png")

