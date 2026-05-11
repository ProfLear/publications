import marimo

__generated_with = "0.19.9"
app = marimo.App(width="medium")


@app.cell
def _():

    from pathlib import Path 
    from PIL import Image
    import numpy as np

    def get_black_pixel_percentage(image_path):
        # Open the image and convert to RGB to ensure standard color channels
        img = Image.open(image_path).convert('RGB')
        pixels = img.getdata()
    
        total_pixels = len(pixels)
        # Define "black" as (0, 0, 0)
        black_pixels = sum(1 for pixel in pixels if pixel == (0, 0, 0))
    
        percentage = (black_pixels / total_pixels) * 100
        return percentage



    return Path, get_black_pixel_percentage, np


@app.cell
def _(Path, get_black_pixel_percentage, np):
    # Usage

    for pattern in ["flat", "small", "large"]:
        sample_path = Path(f"./{pattern}/week3")
        percent = []
    
        for n in ["a", "b", "c", "d"]:
            path = sample_path/f"{pattern}_week3{n}_segmented.png"
        
            try:
                result = get_black_pixel_percentage(path)
                percent.append(result)
                print(f"for {pattern}_{n}, percent uncovered: {result:.1f}%")
            except FileNotFoundError:
                print("Error: Image file not found.")
            
        #print(percent)
        percent = np.array(percent)
        print(f"average = {np.mean(percent):.1f}, stderr = {np.std(percent)/np.sqrt(len(percent)):.1f}")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
