import glob
from plotly.subplots import make_subplots
import numpy as np

R_max_grid_dict = {}
R_max_grid_dict["name"] = []
R_max_grid_dict["R_max"] = []
R_max_grid_dict["AIC"] = []
R_max_grid_dict["I_scale"] = []
R_max_grid_dict["ell"] = []
R_max_grid_dict["blindspot"] = []
R_max_grid_dict["width_scale"] = []

file_pattern = 'MimsFitting/R_max_grid_search/*.csv'

R_max_files = glob.glob(file_pattern, recursive=True)

for _file in R_max_files:
    R_max_grid_dict["name"].append(_file.split("/")[-1].split(".")[0])
    with open(_file, "r") as _f:
        for _line in _f:
            for _key in R_max_grid_dict:
                if _key in _line.split(",")[0]:
                    R_max_grid_dict[_key].append(float(_line.split(",")[1]))

R_max_plot = make_subplots(rows = len(R_max_grid_dict)-2 , cols =1)
i = 1
for _key in R_max_grid_dict:
    if _key != "R_max" and _key != "name":
        R_max_plot.add_scatter(x = R_max_grid_dict["R_max"], y = R_max_grid_dict[_key], name = _key, mode = "markers", row = i, col = 1)
        i = i+1
        if _key == "blindspot": 
            for _j, _item in enumerate(R_max_grid_dict[_key]):
                if _item < 0.05:
                    print(f"{R_max_grid_dict['name'][_j]}: AIC = {R_max_grid_dict['AIC'][_j]}, ell = {R_max_grid_dict['ell'][_j]}")
                    


R_max_plot.update_xaxes(
    range = [0, max(R_max_grid_dict["R_max"])*1.02],
    )
R_max_plot.show()

print(np.mean(np.array(R_max_grid_dict["ell"])))
print(np.mean(np.array(sorted(R_max_grid_dict["ell"][int(len(R_max_grid_dict["ell"])*0.25):int(len(R_max_grid_dict["ell"])*0.75)]))))
print(np.median(np.array(sorted(R_max_grid_dict["ell"][int(len(R_max_grid_dict["ell"])*0.25):int(len(R_max_grid_dict["ell"])*0.75)]))))

print(np.min(np.array(R_max_grid_dict["I_scale"])))
print(np.min(np.array(R_max_grid_dict["blindspot"])))

