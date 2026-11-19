# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "marimo",
#     "numpy",
#     "pillow",
#     "matplotlib",
# ]
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import io
    import os
    from pathlib import Path
    import marimo as mo
    import matplotlib.colors as mcolors
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont

    return Image, ImageDraw, ImageFont, Path, io, mcolors, mo, np, os


@app.cell
def _(Path, mo):
    mo.md(
        """
        # Interactive HSV Channel Fixing & Luminosity Notebook
        Select an image and use the **single set of HSV sliders** below to choose the fixed values $(H_{\\text{fix}}, S_{\\text{fix}}, V_{\\text{fix}})$ applied whenever a channel is held constant across the 7 transformations:
        1. **S and V are fixed, H varies** — $(H(x,y),\\, S_{\\text{fix}},\\, V_{\\text{fix}})$
        2. **H and V are fixed, S varies** — $(H_{\\text{fix}},\\, S(x,y),\\, V_{\\text{fix}})$
        3. **H and S are fixed, V varies** — $(H_{\\text{fix}},\\, S_{\\text{fix}},\\, V(x,y))$
        4. **S and V vary, H is fixed** — $(H_{\\text{fix}},\\, S(x,y),\\, V(x,y))$
        5. **H and V vary, S is fixed** — $(H(x,y),\\, S_{\\text{fix}},\\, V(x,y))$
        6. **H and S vary, V is fixed** — $(H(x,y),\\, S(x,y),\\, V_{\\text{fix}})$
        7. **Black and white (luminosity)** — Grayscale luminosity ($Y = 0.2126 R + 0.7152 G + 0.0722 B$)
        """
    )
    script_dir = Path(__file__).parent.resolve()
    output_dir = script_dir / "output"
    output_dir.mkdir(exist_ok=True)
    return output_dir, script_dir


@app.cell
def _(Path, mo, script_dir):
    # Discover images in HSVmodifications/ and ../OverleafProject/Figures/
    image_extensions = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".webp"}
    discovered_files = {}

    for p in sorted(script_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in image_extensions:
            discovered_files[f"HSVmodifications/{p.name}"] = str(p)

    figures_dir = script_dir.parent / "OverleafProject" / "Figures"
    if figures_dir.exists():
        for p in sorted(figures_dir.iterdir()):
            if p.is_file() and p.suffix.lower() in image_extensions:
                discovered_files[f"Figures/{p.name}"] = str(p)

    default_key = (
        "HSVmodifications/lavoisier_met_crop.png"
        if "HSVmodifications/lavoisier_met_crop.png" in discovered_files
        else (next(iter(discovered_files.keys())) if discovered_files else None)
    )

    preset_dropdown = mo.ui.dropdown(
        options=discovered_files,
        value=default_key,
        label="Select workspace image:",
        full_width=True,
    )

    uploaded_file = mo.ui.file(
        filetypes=[".png", ".jpg", ".jpeg", ".tif", ".tiff", ".webp"],
        label="Or upload/drop any image file (overrides dropdown when uploaded):",
        multiple=False,
    )

    preview_max_dim = mo.ui.dropdown(
        options={"Fast (600 px)": 600, "Medium (1000 px)": 1000, "Full Resolution": 0},
        value="Fast (600 px)",
        label="Interactive preview resolution (Export is always full resolution):",
    )

    crop_left = mo.ui.slider(
        start=0.0,
        stop=90.0,
        step=0.5,
        value=0.0,
        label="Crop Left (%):",
        show_value=True,
        full_width=True,
    )
    crop_right = mo.ui.slider(
        start=0.0,
        stop=90.0,
        step=0.5,
        value=0.0,
        label="Crop Right (%):",
        show_value=True,
        full_width=True,
    )
    crop_top = mo.ui.slider(
        start=0.0,
        stop=90.0,
        step=0.5,
        value=0.0,
        label="Crop Top (%):",
        show_value=True,
        full_width=True,
    )
    crop_bottom = mo.ui.slider(
        start=0.0,
        stop=90.0,
        step=0.5,
        value=0.0,
        label="Crop Bottom (%):",
        show_value=True,
        full_width=True,
    )

    mo.vstack(
        [
            mo.md("### 1. Select Input Image & Cropping (%)"),
            mo.hstack([preset_dropdown, uploaded_file, preview_max_dim], widths="equal"),
            mo.hstack([crop_left, crop_right, crop_top, crop_bottom], widths="equal"),
        ]
    )
    return (
        crop_bottom,
        crop_left,
        crop_right,
        crop_top,
        default_key,
        discovered_files,
        figures_dir,
        image_extensions,
        preset_dropdown,
        preview_max_dim,
        uploaded_file,
    )


@app.cell
def _(mo):
    h_slider = mo.ui.slider(
        start=0.0,
        stop=1.0,
        step=0.01,
        value=0.08,
        label="Fixed Hue (H):",
        show_value=True,
        full_width=True,
    )
    s_slider = mo.ui.slider(
        start=0.0,
        stop=1.0,
        step=0.01,
        value=1.0,
        label="Fixed Saturation (S):",
        show_value=True,
        full_width=True,
    )
    v_slider = mo.ui.slider(
        start=0.0,
        stop=1.0,
        step=0.01,
        value=1.0,
        label="Fixed Value (V):",
        show_value=True,
        full_width=True,
    )
    reverse_s = mo.ui.checkbox(
        value=False,
        label="Reverse S (1 − S)",
    )
    reverse_v = mo.ui.checkbox(
        value=False,
        label="Reverse V (1 − V)",
    )
    return h_slider, reverse_s, reverse_v, s_slider, v_slider


@app.cell
def _(h_slider, mcolors, mo, np, reverse_s, reverse_v, s_slider, v_slider):
    h_val = float(h_slider.value)
    s_val = float(s_slider.value)
    v_val = float(v_slider.value)
    rev_s = bool(reverse_s.value)
    rev_v = bool(reverse_v.value)

    swatch_rgb = mcolors.hsv_to_rgb([h_val, s_val, v_val])
    swatch_hex = mcolors.to_hex(swatch_rgb)
    deg = int(round(h_val * 360))
    s_pct = int(round(s_val * 100))
    v_pct = int(round(v_val * 100))

    rev_notes = []
    if rev_s:
        rev_notes.append("S reversed (1 − S)")
    if rev_v:
        rev_notes.append("V reversed (1 − V)")
    rev_str = f" &nbsp;|&nbsp; <strong>{', '.join(rev_notes)}</strong>" if rev_notes else ""

    swatch_html = mo.Html(
        f"""
        <div style="display:flex; align-items:center; gap:16px; padding:10px 14px; border:1px solid #ccc; border-radius:8px; background:#fafafa;">
            <div style="width:56px; height:56px; border-radius:6px; background:{swatch_hex}; border:1px solid #333; flex-shrink:0;"></div>
            <div style="font-size:14px; line-height:1.4;">
                <strong>Fixed Reference Color</strong><br/>
                <code>H = {h_val:.2f} ({deg}°)</code> &nbsp;|&nbsp;
                <code>S = {s_val:.2f} ({s_pct}%)</code> &nbsp;|&nbsp;
                <code>V = {v_val:.2f} ({v_pct}%)</code> &nbsp;|&nbsp;
                <code>Hex: {swatch_hex}</code>{rev_str}
            </div>
        </div>
        """
    )

    mo.vstack(
        [
            mo.md("### 2. Fixed HSV Parameter Sliders & Channel Reversal"),
            mo.hstack(
                [
                    mo.vstack([h_slider]),
                    mo.vstack([s_slider, reverse_s]),
                    mo.vstack([v_slider, reverse_v]),
                ],
                widths="equal",
            ),
            swatch_html,
        ]
    )
    return (
        deg,
        h_val,
        rev_notes,
        rev_s,
        rev_str,
        rev_v,
        s_pct,
        s_val,
        swatch_hex,
        swatch_html,
        swatch_rgb,
        v_pct,
        v_val,
    )


@app.cell
def _(
    Image,
    Path,
    crop_bottom,
    crop_left,
    crop_right,
    crop_top,
    io,
    mcolors,
    np,
    preset_dropdown,
    preview_max_dim,
    uploaded_file,
):
    # Load full-resolution image from upload or dropdown
    if uploaded_file.value:
        raw_bytes = uploaded_file.contents()
        raw_pil = Image.open(io.BytesIO(raw_bytes)).convert("RGB")
        source_name = uploaded_file.name() or "uploaded_image"
    elif preset_dropdown.value:
        img_path = Path(preset_dropdown.value)
        raw_pil = Image.open(img_path).convert("RGB")
        source_name = img_path.name
    else:
        raise ValueError("Please select or upload an image.")

    orig_w, orig_h = raw_pil.size

    # Apply percentage cropping (Left, Right, Top, Bottom)
    cl = float(crop_left.value) / 100.0
    cr = float(crop_right.value) / 100.0
    ct = float(crop_top.value) / 100.0
    cb = float(crop_bottom.value) / 100.0

    x0 = min(int(round(orig_w * cl)), orig_w - 2)
    x1 = max(x0 + 1, int(round(orig_w * (1.0 - cr))))
    y0 = min(int(round(orig_h * ct)), orig_h - 2)
    y1 = max(y0 + 1, int(round(orig_h * (1.0 - cb))))

    full_pil = raw_pil.crop((x0, y0, x1, y1))
    full_w, full_h = full_pil.size

    # Prepare working image for interactive display
    max_dim = preview_max_dim.value
    if max_dim and max(full_w, full_h) > max_dim:
        scale = max_dim / float(max(full_w, full_h))
        work_w = max(1, int(round(full_w * scale)))
        work_h = max(1, int(round(full_h * scale)))
        work_pil = full_pil.resize((work_w, work_h), Image.Resampling.LANCZOS)
    else:
        work_pil = full_pil
        work_w, work_h = full_w, full_h

    def compute_hsv_variants(pil_img, h_fix, s_fix, v_fix, invert_s=False, invert_v=False):
        rgb = np.asarray(pil_img, dtype=np.float32) / 255.0
        hsv = mcolors.rgb_to_hsv(rgb)
        H = hsv[:, :, 0]
        S = 1.0 - hsv[:, :, 1] if invert_s else hsv[:, :, 1]
        V = 1.0 - hsv[:, :, 2] if invert_v else hsv[:, :, 2]

        H_fixed = np.full_like(H, h_fix)
        S_fixed = np.full_like(S, s_fix)
        V_fixed = np.full_like(V, v_fix)

        def to_pil(h_ch, s_ch, v_ch):
            arr_rgb = np.clip(mcolors.hsv_to_rgb(np.dstack([h_ch, s_ch, v_ch])), 0.0, 1.0)
            return Image.fromarray((arr_rgb * 255.0).round().astype(np.uint8))

        # 1. S and V are fixed, H varies
        img_1 = to_pil(H, S_fixed, V_fixed)
        # 2. H and V are fixed, S varies
        img_2 = to_pil(H_fixed, S, V_fixed)
        # 3. H and S are fixed, V varies
        img_3 = to_pil(H_fixed, S_fixed, V)
        # 4. S and V vary, H is fixed
        img_4 = to_pil(H_fixed, S, V)
        # 5. H and V vary, S is fixed
        img_5 = to_pil(H, S_fixed, V)
        # 6. H and S vary, V is fixed
        img_6 = to_pil(H, S, V_fixed)
        # 7. Black and white (luminosity: ITU-R BT.709 relative luminance)
        Y = np.clip(0.2126 * rgb[:, :, 0] + 0.7152 * rgb[:, :, 1] + 0.0722 * rgb[:, :, 2], 0.0, 1.0)
        img_7 = Image.fromarray((np.dstack([Y, Y, Y]) * 255.0).round().astype(np.uint8))

        return {
            "00_original": pil_img,
            "01_SV_fixed_H_varies": img_1,
            "02_HV_fixed_S_varies": img_2,
            "03_HS_fixed_V_varies": img_3,
            "04_SV_vary_H_fixed": img_4,
            "05_HV_vary_S_fixed": img_5,
            "06_HS_vary_V_fixed": img_6,
            "07_black_and_white_luminosity": img_7,
        }

    return (
        compute_hsv_variants,
        full_h,
        full_pil,
        full_w,
        max_dim,
        source_name,
        work_h,
        work_pil,
        work_w,
    )


@app.cell
def _(
    compute_hsv_variants,
    deg,
    full_h,
    full_w,
    h_val,
    mo,
    rev_s,
    rev_v,
    s_pct,
    s_val,
    source_name,
    v_pct,
    v_val,
    work_pil,
):
    variants_preview = compute_hsv_variants(
        work_pil, h_val, s_val, v_val, invert_s=rev_s, invert_v=rev_v
    )

    s_tag = "1−S varies" if rev_s else "S varies"
    v_tag = "1−V varies" if rev_v else "V varies"

    labels_map = [
        ("00_original", "Original Image", f"{source_name} ({full_w}×{full_h})"),
        ("07_black_and_white_luminosity", "7. Black & White (Luminosity)", "Y = 0.2126R + 0.7152G + 0.0722B"),
        ("01_SV_fixed_H_varies", "1. S & V fixed, H varies", f"H varies | S={s_pct}%, V={v_pct}%"),
        ("02_HV_fixed_S_varies", "2. H & V fixed, S varies", f"{s_tag} | H={deg}°, V={v_pct}%"),
        ("03_HS_fixed_V_varies", "3. H & S fixed, V varies", f"{v_tag} | H={deg}°, S={s_pct}%"),
        ("04_SV_vary_H_fixed", "4. S & V vary, H fixed", f"{s_tag}, {v_tag} | H={deg}°"),
        ("05_HV_vary_S_fixed", "5. H & V vary, S fixed", f"H, {v_tag} | S={s_pct}%"),
        ("06_HS_vary_V_fixed", "6. H & S vary, V fixed", f"H, {s_tag} | V={v_pct}%"),
    ]

    def make_card(key, title, subtitle):
        return mo.vstack(
            [
                mo.md(f"**{title}**  \n<small>`{subtitle}`</small>"),
                mo.image(variants_preview[key], width=320),
            ],
            align="center",
        )

    row_0 = mo.hstack(
        [make_card(*labels_map[0]), make_card(*labels_map[1])],
        justify="start",
        gap=2,
    )
    row_1 = mo.hstack(
        [make_card(*labels_map[2]), make_card(*labels_map[3]), make_card(*labels_map[4])],
        justify="start",
        gap=2,
    )
    row_2 = mo.hstack(
        [make_card(*labels_map[5]), make_card(*labels_map[6]), make_card(*labels_map[7])],
        justify="start",
        gap=2,
    )

    mo.vstack(
        [
            mo.md("### 3. Live Preview of All 7 Recolored Images (+ Original)"),
            mo.md("#### Original & Black-and-White Luminosity"),
            row_0,
            mo.md("#### Two Channels Fixed, One Channel Varies (1–3)"),
            row_1,
            mo.md("#### One Channel Fixed, Two Channels Vary (4–6)"),
            row_2,
        ]
    )
    return labels_map, make_card, row_0, row_1, row_2, s_tag, v_tag, variants_preview


@app.cell
def _(mo, output_dir):
    save_button = mo.ui.run_button(
        label=f"Save Full-Resolution Images to {output_dir.name}/",
        kind="success",
    )
    mo.vstack(
        [
            mo.md("### 4. Export Full-Resolution Images"),
            save_button,
        ]
    )
    return (save_button,)


@app.cell
def _(
    Image,
    ImageDraw,
    ImageFont,
     Path,
    compute_hsv_variants,
    full_pil,
    h_val,
    mo,
    output_dir,
    rev_s,
    rev_v,
    s_val,
    save_button,
    v_val,
):
    # Always save or update when save_button is clicked (or on initial run if output doesn't have the 7 files)
    should_save = save_button.value or not (output_dir / "01_SV_fixed_H_varies.png").exists()

    saved_paths = []
    if should_save:
        full_variants = compute_hsv_variants(
            full_pil, h_val, s_val, v_val, invert_s=rev_s, invert_v=rev_v
        )
        for k, im in full_variants.items():
            out_file = output_dir / f"{k}.png"
            im.save(out_file, dpi=(300, 300))
            saved_paths.append(str(out_file))

        # Also build a clean composite sheet
        fw, fh = full_pil.size
        target_w = 480
        target_h = int(round(target_w * (fh / fw)))

        def prep_panel(im):
            p = im.resize((target_w, target_h), Image.Resampling.LANCZOS)
            d = ImageDraw.Draw(p)
            d.rectangle([0, 0, target_w - 1, target_h - 1], outline=(50, 50, 50), width=1)
            return p

        header_h = 50
        pad_x = 24
        pad_y = 45
        margin_left = 60
        margin_right = 25
        margin_bottom = 25

        # 4 columns x 2 rows (Left column: Original on top, B&W Luminosity on bottom; Right 3 columns: 1-3 on top, 4-6 on bottom)
        total_w = margin_left + 4 * target_w + 3 * pad_x + margin_right
        total_h = header_h + 2 * target_h + pad_y + margin_bottom

        comp = Image.new("RGB", (total_w, total_h), color=(255, 255, 255))
        draw = ImageDraw.Draw(comp)
        try:
            font_header = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 19)
        except Exception:
            font_header = ImageFont.load_default()

        grid_items = [
            # Top row
            (0, 0, "(a) original", full_variants["00_original"]),
            (1, 0, "(b) S & V fixed, H varies", full_variants["01_SV_fixed_H_varies"]),
            (2, 0, "(c) H & V fixed, S varies", full_variants["02_HV_fixed_S_varies"]),
            (3, 0, "(d) H & S fixed, V varies", full_variants["03_HS_fixed_V_varies"]),
            # Bottom row
            (0, 1, "(e) black & white (luminosity)", full_variants["07_black_and_white_luminosity"]),
            (1, 1, "(f) S & V vary, H fixed", full_variants["04_SV_vary_H_fixed"]),
            (2, 1, "(g) H & V vary, S fixed", full_variants["05_HV_vary_S_fixed"]),
            (3, 1, "(h) H & S vary, V fixed", full_variants["06_HS_vary_V_fixed"]),
        ]

        for col_idx, row_idx, title_text, panel_im in grid_items:
            px = margin_left + col_idx * (target_w + pad_x)
            py = header_h + row_idx * (target_h + pad_y)
            comp.paste(prep_panel(panel_im), (px, py))
            draw.text((px, py - 28), title_text, fill=(0, 0, 0), font=font_header)

        comp_path = output_dir / "composite_hsv_fixed.png"
        comp.save(comp_path, dpi=(300, 300))
        saved_paths.append(str(comp_path))

    status_md = (
        mo.md(
            f"**Saved {len(saved_paths)} full-resolution files** to `{output_dir}` (with $H={h_val:.2f}, S={s_val:.2f}, V={v_val:.2f}$, `reverse_S={rev_s}`, `reverse_V={rev_v}`):\n"
            + "\n".join(f"- `{Path(p).name}`" for p in saved_paths)
        )
        if saved_paths
        else mo.md("Click the button above to export full-resolution images with the current slider values.")
    )
    status_md
    return comp_path, full_variants, saved_paths, should_save, status_md


if __name__ == "__main__":
    app.run()
