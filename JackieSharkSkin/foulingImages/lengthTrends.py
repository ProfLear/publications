import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from plotly.subplots import make_subplots
    import numpy as np
    from PIL import Image
    from codechembook import quickPlots as qp
    from pathlib import Path

    return Image, Path, make_subplots, np, qp


@app.cell
def _(Image, np):
    def get_class_coverage(image_path):
        img = Image.open(image_path).convert("RGB")
        pixels = img.getdata()

        total_pixels = len(pixels)
        clear_pixels = 0
        thick_pixels = 0
        for pixel in pixels:
            if pixel[0] > 200 and pixel[1] < 10 and pixel[2] < 10:
                clear_pixels +=1
            elif pixel[0] > 100 and pixel[0] < 200 and pixel[1] > 100 and pixel[1] < 200 and pixel[2] > 200: # this is purple
                thick_pixels += 1


        thin_pixels = total_pixels - clear_pixels - thick_pixels

        return np.array([clear_pixels, thin_pixels, thick_pixels])/total_pixels

    return (get_class_coverage,)


@app.cell
def _(Path, get_class_coverage):
    patterns = ["flat", "small", "large"]
    samples = ["a", "b", "c", "d"]
    colors = [(0, 54, 79), (79, 0, 45), (12, 77, 0)]
    all_coverage = {} # Renamed to avoid collision
    root = "retraining_fouling" # blank for mine or /jackie_weka
    for _pattern in patterns:
        all_coverage[_pattern] = {}
        for _week in range(1, 13):
            folder = Path(f"./{root}/{_pattern}/week{_week}")
            for _sample in samples:
                image_path = folder / f"{_pattern}_{_week}{_sample}_segmented.jpg"

                try:
                    # Use a unique name for the function output
                    temp_coverage = get_class_coverage(image_path) 
                    if _week not in all_coverage[_pattern]:
                        all_coverage[_pattern][_week] = {}
                    all_coverage[_pattern][_week][_sample] = temp_coverage
                    print(f"{image_path} read!")
                except Exception as e:
                    print(f"{image_path} not found! Error: {e}")
    return


@app.cell
def _(Image, np):
    def get_linear_trend(image_path, transposed = False):
        img = Image.open(image_path).convert('RGB')
        data = np.array(img)
        if transposed:
            data = np.transpose(data, (1, 0, 2)) #transpose the x and y directions, preserving the triplet rgb
        #print(data)
        rows = []
        fractions_clear = []
        for i, row in enumerate(data):
            clear_pixels = 0
            for pixel in row:
                #print(pixel)
                if pixel[0] > 200 and pixel[1] < 10 and pixel[2] < 10:
                    clear_pixels +=1
            rows.append(i)
            fractions_clear.append(clear_pixels/len(row))

        return rows, fractions_clear

    def get_radial_trend(image_path, bin_size = 1):
        img = Image.open(image_path).convert('RGB')
        data = np.array(img)
        h, w, _ = data.shape
    
        # 1. Define the center
        center_y, center_x = h // 2, w // 2
    
        # 2. Create a grid of coordinates
        y, x = np.indices((h, w))
        r = np.sqrt((x - center_x)**2 + (y - center_y)**2).astype(int)

        # 3. Create a mask for "clear" pixels (Red > 200, G < 10, B < 10)
        # This is much faster than a nested for-loop
        clear_mask = (data[:,:,0] > 200) & (data[:,:,1] < 10) & (data[:,:,2] < 10)

        # 4. Analyze by radius
        max_radius = min(center_x, center_y)
        radii = np.arange(0, max_radius, bin_size)
        fractions_clear = []

        for radius in radii:
            # Find all pixels that fall into this specific radial "ring"
            mask_at_r = (r == radius)
            pixels_at_r = clear_mask[mask_at_r]
        
            if len(pixels_at_r) > 0:
                fraction = np.sum(pixels_at_r) / len(pixels_at_r)
            else:
                fraction = 0
            
            fractions_clear.append(fraction)

        return radii.tolist(), fractions_clear
    
    def get_stats_truncated(arrays):
        # 1. Find the shortest length
        min_len = min(len(a) for a in arrays)
    
        # 2. Truncate all arrays
        truncated = [a[:min_len] for a in arrays]
    
        # 3. Calculate mean and std along axis 0
        mean_array = np.mean(truncated, axis=0)
        std_array = np.std(truncated, axis=0)
    
        return mean_array, std_array

    return get_linear_trend, get_radial_trend, get_stats_truncated


@app.cell
def _(get_radial_trend, get_stats_truncated, make_subplots):


    radial_plot = make_subplots(rows = 3, cols = 1, subplot_titles=["flat", "small", "large"])
    for _row, _pattern in enumerate(["flat", "small", "large"], start = 1):
        collected = []
        for _sample in ["a", "b", "c", "d"]:
            r, frac_clear = get_radial_trend(f"/Users/benjaminlear/GitHub/publications/JackieSharkSkin/retraining_fouling/{_pattern}/washed/{_pattern}_washed{_sample}_segmented.jpg", bin_size = 1)
            radial_plot.add_scatter(x = r, y = frac_clear,
                                    line = dict(color = "lightgrey"),
                                    showlegend=False,
                                    row  = _row, col = 1)
            collected.append(frac_clear)
        mean, std = get_stats_truncated(collected)
    
        radial_plot.add_scatter(x = r, y = mean+std/2, 
                                    line = dict(color = "black", width = 0),
                                    showlegend=False,
                                    row  = _row, col = 1)
        radial_plot.add_scatter(x = r, y = mean - std/2, 
                                    line = dict(color = "black", width = 0),
                                    fill = "tonexty",
                                    showlegend=False,
                                    row  = _row, col = 1)
        radial_plot.add_scatter(x = r, y = mean, 
                                    line = dict(color = "black"),
                                    showlegend=False,
                                    row  = _row, col = 1)

    radial_plot.update_xaxes(title = "distance from center / pixels", row = 3, col = 1)
    radial_plot.update_yaxes(title = "fraction clear", range = [0, 1])
    radial_plot.update_layout(template = "simple_white", height = 700)

    radial_plot.show()
    return


@app.cell
def _(get_linear_trend, make_subplots, np):
    test_plot = make_subplots()

    x, y = get_linear_trend(r"/Users/benjaminlear/GitHub/publications/JackieSharkSkin/retraining_fouling/large/washed/large_washedb_segmented.jpg", transposed=False)

    shuffle = False
    if shuffle:
        rng = np.random.default_rng(seed=42)
        rng.shuffle(x)


    test_plot.add_scatter(x = x, y = y, mode = "markers")
    test_plot.add_scatter(x = [0, len(x)], y = [np.mean(y)]*2)
    test_plot.add_scatter(x = [0, len(x)], y = [np.mean(y)+np.std(y)]*2)

    test_plot.update_layout(template = "simple_white")
    test_plot.show()
    return (y,)


@app.cell
def _(np, qp, y):
    yf = np.fft.rfft(y)
    xf = np.fft.rfftfreq(len(y), d=1)

    # 3. Get the Magnitude (the "strength" of each frequency)
    amplitude = np.abs(yf)

    qp.quickScatter(x = xf, y = amplitude)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
