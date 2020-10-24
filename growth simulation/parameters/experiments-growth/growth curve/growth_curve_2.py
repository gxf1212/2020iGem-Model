import numpy as np
import xlrd
import os
import time
from utils import *

# %% organize files
path = r'../growth/'
files_path = list(os.listdir(path))  # list all paths
files_path.remove('empty.xlsx')  # not using
files_time = []  # using start time
for fpath in files_path:
    data = xlrd.open_workbook(path + fpath).sheet_by_index(0)
    ftime = data.cell_value(21,1)
    # ftime = time.ctime(os.stat(path+fpath).st_mtime)[8:19]  # find file modifying time
    files_time.append(t2s(ftime))
z = zip(files_path, files_time)
zipped = sorted(z, key=lambda z: z[1])  # sort by time
files_path, files_time = zip(*zipped)  # unzip

# %% read data
time = (files_time-np.min(files_time))/3600  # convert to hour
value = []
for file in files_path:
    data = xlrd.open_workbook(path+file).sheet_by_index(0)
    v = []
    for i in [26, 27, 28]:  # row number
        v.append(data.row_values(i-1, 1, 9))
    v_ave = np.mean(np.array(v), axis=0)  # mean for each file
    value.append(v_ave)
value = np.array(value)
# check for abnormal value
# plot1(value, time)
# linear interpolation for index 13,15,23
for i in [13,15,23]:
    for j in range(value.shape[1]):
        value[i,j]=(value[i-1,j]+value[i+1,j])/2
plot1(value, time)
# each column minus empty control group
value = np.array([value[:,i]-value[:,7] for i in range(7)])

#%% plot
plt.figure()
for i in range(value.shape[0]):
    plt.plot(time, value[i,:])
plt.legend(["1","2","3","4","5","6","7"])
plt.show()

# %% save
np.savetxt('time.csv', time, delimiter=',')
np.savetxt('growth.csv', value, delimiter=',')

