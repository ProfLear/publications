import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def import_and_downsample(xyz_file, sampling = 1):
    with open (xyz_file, "r") as xyz:
        # this block gets the information on the image width and height
        n_headers = 0
        n_data_rows = 0
        n_data_cols = 0
        zs = []
        for i, line in enumerate(xyz, start = 1):
            line = line.strip("\n")  # let us remove newline characters
            if line[0] == "#": # we are in a header
                n_headers = n_headers + 1
                if "Width" in line:
                    total_width = float(line.split(" ")[-2])
                    #print(line.split(" "))
                    if line.split(" ")[-1] == "mm":
                        total_width = total_width * 1E-3 # convert from mm to meters
                    else:
                        print("width unit is not recognized")
                if "Height" in line:
                    total_height = float(line.split(" ")[-2])
                    if line.split(" ")[-1] == "mm":
                        total_height = total_height * 1E-3 # convert from mm to meters
                    else:
                        print("height unit is not recognized")
                if "Value units" in line:
                    if line.split(" ")[-1] == "m":
                        z_scale = 1 # what we should multiply by to get to meters

            else: # then we are not in a header. This is data. 
                n_data_rows = n_data_rows + 1
                if (i - n_headers)%sampling == 0:# this is a row we will sample
                    temp_row = []
                    for j, value in enumerate(line.split("\t")):
                        if j%sampling == 0: # this is a point to sample
                            temp_row.append(float(value)*z_scale)
                            final_sampled_col = j
                    if n_data_cols == 0: # we haven't yet recorded the number of data points
                        n_data_cols = j # so recordit
                    zs.append(np.array(temp_row))
                    final_sampled_row = i

    zs = np.array(zs)
    return zs, [(total_width*final_sampled_col/n_data_cols)/len(zs[0]), (total_height*final_sampled_row/n_data_rows)/len(zs)] # also return x_scale and y_scale in xy_scales


def plot_single_areal(areal_zs,  xy_scales = [1,1], z_mag = "auto"):
    fig = make_subplots(specs=[[{'type': 'surface',}]])
    fig.add_trace(go.Surface(
        z=areal_zs - np.mean(areal_zs),
        x=np.linspace(0, len(areal_zs), len(areal_zs))*xy_scales[0],
        y=np.linspace(0, len(areal_zs[0]), len(areal_zs[0]))*xy_scales[1],
        colorscale = "greys_r",
    ),
        )
    fig.update_layout(
            autosize=True, 
            margin=dict(l=65, r=50, b=65, t=90),
                    )
    if z_mag != "auto":
        fig.update_layout(
                          margin=dict(l=65, r=50, b=65, t=90),
                          )
    return fig


def make_profile_image(zdata, cutoff = 100, colorscale = "greys_r"):
    profile_image = make_subplots()
    profile_image.add_heatmap(z = zdata, 
        colorscale=colorscale, 
        zmin = np.percentile(zdata, 100-cutoff), 
        zmax= np.percentile(zdata, cutoff),
        colorbar = dict(len = 1.06),
        )
    profile_image.update_xaxes(ticks = "", 
        showticklabels = False,
        showline = True, linewidth = 2, linecolor = "black", mirror = True)
    profile_image.update_yaxes(ticks = "", showticklabels = False,
        showline = True, linewidth = 2, linecolor = "black", mirror = True)
    profile_image.update_layout(width = 500, height = 500)
    return profile_image


def map_paired_heights(map, pairs, low_color = "cyan", high_color = "red", line_color = "white"):
    xs0, ys0, zs0, xs1, ys1, zs1 = [], [], [], [], [], []
    for pair in pairs:
        xs0.append(pair[0][0])
        ys0.append(pair[0][1])
        zs0.append(pair[0][2])

        xs1.append(pair[1][0])
        ys1.append(pair[1][1])
        zs1.append(pair[1][2])

    xs0 = np.array(xs0)
    ys0 = np.array(ys0)
    zs0 = np.array(zs0)

    xs1 = np.array(xs1)
    ys1 = np.array(ys1)
    zs1 = np.array(zs1)


    xs_to_plot = []
    ys_to_plot = []
    colors = []

    for x0, x1, y0, y1 in (zip(xs0, xs1, ys0, ys1)):
        xs_to_plot.append(x0)
        xs_to_plot.append(x1)
        xs_to_plot.append("NA")

        ys_to_plot.append(y0)
        ys_to_plot.append(y1)
        ys_to_plot.append("NA")

        colors.append(low_color)
        colors.append(high_color)
        colors.append("orange")

    map.add_scatter(x = xs_to_plot, y = ys_to_plot, mode = "lines+markers", 
        line = dict(width = 2, color  = line_color),
        marker = dict(size = 6, color = colors, line = dict(width = 2, color = line_color)),
        showlegend = False,
        )


    map.update_layout(showlegend=False)
    #print(zs0)
    #print(zs1)
    print(f"mean height difference = {np.abs(np.mean(zs0-zs1))}, std = {np.std(zs0-zs1)}")

    return map