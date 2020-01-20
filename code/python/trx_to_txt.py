### trx_to_txt.py

import glob, os
import scipy.io as spio

dirpath = "/home/alastair/Downloads/For_Jane/For_Jane/JB"

pathsearch = os.path.join(dirpath,"**","*.mat")
for f in glob.glob(pathsearch,recursive=True):
    print(f)


f

data = spio.loadmat(f, squeeze_me=True, struct_as_record=False)

idx = 5

et = data["trx"][idx].t
run = data["trx"][idx].run

plt.plot(et,run)
plt.show()


dir(data["trx"][idx])