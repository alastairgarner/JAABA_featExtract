
### Dependencies
import scipy.io as spio
import matplotlib.pyplot as plt


### Functions
# def loadmat(filename):
#     '''
#     this function should be called instead of direct spio.loadmat
#     as it cures the problem of not properly recovering python dictionaries
#     from mat files. It calls the function check keys to cure all entries
#     which are still mat-objects
#     '''
#     data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
#     return _check_keys(data)

# def _check_keys(dict):
#     '''
#     checks if entries in dictionary are mat-objects. If yes
#     todict is called to change them to nested dictionaries
#     '''
#     for key in dict:
#         if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
#             dict[key] = _todict(dict[key])
#     return dict        

# def _todict(matobj):
#     '''
#     A recursive function which constructs from matobjects nested dictionaries
#     '''
#     dict = {}
#     for strg in matobj._fieldnames:
#         elem = matobj.__dict__[strg]
#         if isinstance(elem, spio.matlab.mio5_params.mat_struct):
#             dict[strg] = _todict(elem)
#         else:
#             dict[strg] = elem
#     return dict


### Scripting

filename = "data_processed/20190423_112850@GMR_SS01817n@UAS_Chrimson_attp18_72F11@r_LED05_30s2x15s30s#n#n#n.mat"
filename = "./data/jb_results/t93/attp2@UAS_Chrimson_attp18_69F06/r_LED30_45s2x30s30s#n#n#n@100/20180206_155058/trx.mat"

data = loadmat(filename)
# data = spio.loadmat(filename, squeeze_me=True, struct_as_record=False,variable_names=["temp"])
data = spio.loadmat(filename, squeeze_me=True, struct_as_record=False)

type(data)

x = data["temp"].data_chore[1].et
y = data["temp"].data_chore[1].curve

fig,ax = plt.subplots()
ax.plot(x,y)
plt.show()


idx = 5

et = data["trx"][idx].t
run = data["trx"][idx].run

plt.plot(et,run)
plt.show()


dir(data["trx"][idx])