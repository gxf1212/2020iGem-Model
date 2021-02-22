import numpy as np
import xlrd
import os
import pandas as pd
import csv
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
value_all = []  # all data, before average
for file in files_path:
    data = xlrd.open_workbook(path+file).sheet_by_index(0)
    v = []
    for i in [26, 27, 28]:  # row number
        v.append(data.row_values(i-1, 1, 9))
    v_ave = np.mean(np.array(v), axis=0)  # mean for each file
    value.append(v_ave)
    value_all.append(np.array(v))
value = np.array(value)
value_all = np.array(value_all)

#%% write data
# workbook = xlrd.open_workbook(r'./data.xlsx')
#     writer.writerow(value_all[:,0,:])

for i in range(8):
    data = pd.DataFrame(value_all[:,:,i])
    data.to_csv('data'+str(i+1)+'.csv',header = False, index = False)

# with open('data.csv','w') as csvfile:
#     writer = csv.writer(csvfile)
#     for i in range(8):
#         for j in range(3):
#             writer.writerow(value_all[:, :, i])
    # np.savetxt('growth.csv', value, delimiter=',')

#%% check for abnormal value
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

