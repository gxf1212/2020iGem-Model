import matplotlib.pyplot as plt
import numpy as np

# %% input data
t = [0, 2, 4, 6, 8, 10, 12, 13, 14, 15, 16 + 7 / 60, 16 + 15 / 60, 17, 18, 20, 22]  # hour
l1 = [0.362, 0.368, 0.395, 0.449, 0.443, 0.435, 0.425, 0.456, 0.439, 0.426, 0.423]
l2 = [0.362, 0.368, 0.395, 0.487, 0.443, 0.435, 0.425, 0.456, 0.439, 0.433, 0.425]
od1 = [0.007, 0.069, 0.021, 0.397, 1.292] + [i * 8 for i in l1]
od2 = [0.007, -0.006, 0.012, 0.406, 1.332] + [i * 8 for i in l2]
# 0.069/-0.006
# 0.834, 0.867 * 2


# %% preprocess data
data = np.vstack((np.array([t]), np.array([od1])))
idx = []

plt.plot(t, od1)# %% fit the overall growth curve

plt.xlabel('time (hours)')
plt.ylabel('OD600')
plt.show()

