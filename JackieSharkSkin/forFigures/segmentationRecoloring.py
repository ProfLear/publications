import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from PIL import Image
    import os

    return Image, os


@app.cell
def _(Image, os):
    def recolor_clear_pixels(input_image_path, output_folder, old_color = (255, 0, 0), new_color=(0, 255, 0), covered_color = (255, 255, 255), threshold=25):
        from pathlib import Path
        """
        Opens a JPG, replaces 'red' pixels with a new color, and saves the result.

        :param input_path: Path to the original JPG.
        :param output_path: Where to save the new image.
        :param new_color: RGB tuple for the replacement color (default is green).
        :param threshold: How 'red' a pixel needs to be to get swapped.
        """
        image_path = Path(input_image_path)
        output_path = Path(output_folder)
        # Open the image and ensure it's in RGB mode
        img = Image.open(image_path).convert("RGB")
        pixels = img.load()
        width, height = img.size

        old_r, old_g, old_b = old_color

        for y in range(height):
            for x in range(width):
                r, g, b = pixels[x, y]

                # Logic: Is Red significantly higher than Green and Blue?
                if r > old_r-threshold and r < old_r + threshold and  g > old_g - threshold and g < old_g + threshold and b > old_b - threshold and b < old_b + threshold:
                    pixels[x, y] = new_color

                else:
                    pixels[x, y] = covered_color

        # Ensure the output directory exists
        os.makedirs(os.path.dirname(output_folder), exist_ok=True)

        # Save the result
        img.save(output_path/image_path.name)

        print(f"Success! Image saved to {output_path}")

    return (recolor_clear_pixels,)


@app.cell
def _():
    # 105% lightness
    flat_color =  (7, 115, 177) #(2, 89, 127)
    small_color = (231,159, 39) # (255, 118, 199)
    large_color = (162, 55, 104)# (255, 234, 52)

    covered_color = (200, 255, 200)
    return covered_color, flat_color, large_color, small_color


@app.cell
def _(
    covered_color,
    flat_color,
    large_color,
    recolor_clear_pixels,
    small_color,
):
    recolor_clear_pixels(
        r"/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publications/JackieSharkSkin/foulingImages/flat/week2/flat_2a_segmented.JPG",
        r"C./",
        new_color = flat_color,  
        covered_color=covered_color,
    )

    recolor_clear_pixels(
        r"/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publications/JackieSharkSkin/foulingImages/large/week2/large_2a_segmented.JPG",
        r"C./",
        new_color = large_color,  
        covered_color=covered_color,
    )

    recolor_clear_pixels(
        r"/Users/benjaminlear/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/publications/JackieSharkSkin/foulingImages/small/week2/small_2a_segmented.JPG",
        r"C./",
        new_color = small_color,
        covered_color=covered_color,
    )
    return


@app.cell
def _(flat_color, large_color, recolor_clear_pixels, small_color):
    patterns = ["flat", "small", "large"]
    samples = ["a", "b", "c", "d"]
    colors = [flat_color, small_color, large_color]
    for _pattern, _color in zip(patterns, colors):
        for _week in range(13):
            for _sample in samples:
                try:
                    recolor_clear_pixels(f"../foulingImages/{_pattern}/week{_week}/{_pattern}_{_week}{_sample}_segmented.JPG", 
                        f"./segmented_{_pattern}",             
                                old_color = (255, 0, 0), 
                                new_color=_color, 
                                covered_color = (225, 225, 225), 
                                threshold=25)

                except:
                    print(f"../foulingImages/{_pattern}/week{_week}/{_pattern}_{_week}{_sample}_segmented.JPG not found!")
    f"../foulingImages/{_pattern}/week{_week}/{_pattern}_{_week}{_sample}_segmented.JPG",
    return


if __name__ == "__main__":
    app.run()
