import os
import glob
import random
import shutil
import re
import rasterio
from rasterio.errors import RasterioIOError
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# --- Configuration ---
# Specific ROI to process
TARGET_ROI_PATH = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_lst_16days/LeLoi_KienXuong_ThaiBinh"
OUTPUT_VIS_SUBFOLDER = "band_time_series_plot" # Will be created inside TARGET_ROI_PATH/visualization/

# Folder names within the TARGET_ROI_PATH
LST_FOLDER_NAME = "lst"
ERA5_FOLDER_NAME = "era5"
NDVI_FOLDER_NAME = "ndvi_infer"  # Assumption: NDVI .tif files are here

# File patterns - **ADJUST NDVI_PATTERN IF NEEDED**
LST_PATTERN = "lst_*.tif"
ERA5_PATTERN = "era5_*.tif"
NDVI_PATTERN = "ndvi_*.tif"  # ASSUMPTION: e.g., ndvi_2023-01-01.tif. Or could be just "*.tif" if dates are different.

BANDS_CONFIG = {
    "NDVI": {"folder": NDVI_FOLDER_NAME, "pattern": NDVI_PATTERN, "band_index": 1},
    "LST": {"folder": LST_FOLDER_NAME, "pattern": LST_PATTERN, "band_index": 1},
    "ERA5_T2M": {"folder": ERA5_FOLDER_NAME, "pattern": ERA5_PATTERN, "band_index": 1},
    "ERA5_SKT": {"folder": ERA5_FOLDER_NAME, "pattern": ERA5_PATTERN, "band_index": 2},
}
# Order of bands in the plot (adjust if NDVI is not first or if other bands are added)
BANDS_ORDER = ["NDVI", "LST", "ERA5_T2M", "ERA5_SKT"]

NUM_IMAGES_PER_BAND = 20 # Number of random images (columns) per band

# --- Helper Functions ---
def extract_date_from_prepared_filename(filename_str):
    """Extracts YYYY-MM-DD from filenames like 'prefix_YYYY-MM-DD.tif' or 'YYYY-MM-DD.tif'"""
    basename = os.path.basename(filename_str)
    # Try to match prefix_YYYY-MM-DD.tif
    match_prefix = re.search(r"^[a-zA-Z0-9]+_(\d{4}-\d{2}-\d{2})\.tif$", basename)
    if match_prefix:
        return match_prefix.group(1)
    
    # Try to match YYYY-MM-DD.tif (if no prefix)
    match_no_prefix = re.search(r"^(\d{4}-\d{2}-\d{2})\.tif$", basename)
    if match_no_prefix:
        return match_no_prefix.group(1)
        
    return os.path.splitext(basename)[0][:10] # Fallback

def get_band_files(roi_base_path, band_config):
    """Gets all tif files for a specific band configuration in the ROI."""
    folder_path = os.path.join(roi_base_path, band_config["folder"])
    if not os.path.isdir(folder_path):
        print(f"  Info: Folder {folder_path} not found for band.")
        return []
    
    search_pattern = os.path.join(folder_path, band_config["pattern"])
    files = glob.glob(search_pattern)
    if not files:
        print(f"  Info: No files found in {folder_path} matching pattern {band_config['pattern']}")
    return sorted(files)

def visualize_single_roi(roi_path, vis_output_dir):
    """Creates and saves a plot for the specified ROI."""
    roi_name = os.path.basename(roi_path)
    print(f"Visualizing ROI: {roi_name}")

    fig, axes = plt.subplots(len(BANDS_ORDER), NUM_IMAGES_PER_BAND, 
                             figsize=(NUM_IMAGES_PER_BAND * 3, len(BANDS_ORDER) * 3.5), # Slightly more height for colorbars
                             squeeze=False)
    
    plt.suptitle(f"Data Visualization for ROI: {roi_name}", fontsize=16, y=1.0)

    for i, band_name_key in enumerate(BANDS_ORDER):
        current_band_config = BANDS_CONFIG[band_name_key]
        band_files = get_band_files(roi_path, current_band_config)
        
        axes[i, 0].set_ylabel(band_name_key, fontsize=10, rotation=0, labelpad=40, ha='right')

        if not band_files:
            for j_ax in range(NUM_IMAGES_PER_BAND):
                axes[i, j_ax].text(0.5, 0.5, 'No Data', ha='center', va='center', transform=axes[i, j_ax].transAxes)
                axes[i, j_ax].axis('off')
            continue

        if len(band_files) > NUM_IMAGES_PER_BAND:
            selected_files = random.sample(band_files, NUM_IMAGES_PER_BAND)
        else:
            selected_files = band_files
        selected_files.sort()

        for j, file_path in enumerate(selected_files):
            ax = axes[i, j]
            ax.axis('off')
            try:
                with rasterio.open(file_path) as src:
                    data = src.read(current_band_config["band_index"])
                    if src.nodata is not None:
                        nodata_val = np.array(src.nodata, dtype=data.dtype).item()
                        data = np.where(data == nodata_val, np.nan, data.astype(float))
                    else:
                        data = data.astype(float)

                    if np.all(np.isnan(data)):
                        ax.text(0.5, 0.5, 'All NoData', ha='center', va='center', transform=ax.transAxes)
                    else:
                        vmin, vmax = np.nanpercentile(data, [2, 98])
                        if vmin == vmax: 
                            vmin = np.nanmin(data) - 0.1 if np.nanmin(data) is not None else -0.1
                            vmax = np.nanmax(data) + 0.1 if np.nanmax(data) is not None else 0.1
                            if vmin == vmax : vmax += 0.1 # Ensure range if still equal
                        
                        im = ax.imshow(data, cmap='viridis', vmin=vmin, vmax=vmax)
                        plt.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
                
                file_date = extract_date_from_prepared_filename(file_path)
                ax.set_title(file_date if file_date else os.path.basename(file_path)[:10], fontsize=8)
            
            except IndexError:
                ax.text(0.5, 0.5, f'Band {current_band_config["band_index"]}\nNot Found', ha='center', va='center', transform=ax.transAxes, fontsize=8)
            except RasterioIOError:
                ax.text(0.5, 0.5, 'Read Error', ha='center', va='center', transform=ax.transAxes)
            except Exception as e:
                ax.text(0.5, 0.5, 'Plot Error', ha='center', va='center', transform=ax.transAxes)
                print(f"    Error plotting {file_path} for band {band_name_key}: {e}")
        
        # Fill remaining axes if fewer than NUM_IMAGES_PER_BAND were plotted
        for j_fill in range(len(selected_files), NUM_IMAGES_PER_BAND):
            axes[i, j_fill].axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.97]) # Adjust rect to give space for suptitle and bottom
    
    output_filename = f"{os.path.basename(roi_path)}_bands_visualization.png"
    output_png_path = os.path.join(vis_output_dir, output_filename)
    try:
        plt.savefig(output_png_path, dpi=150, bbox_inches='tight')
        print(f"  Saved plot to: {output_png_path}")
    except Exception as e:
        print(f"  Error saving plot {output_png_path}: {e}")
    plt.close(fig)

# --- Main Execution ---
def main():
    print(f"Targeting ROI: {TARGET_ROI_PATH}")
    if not os.path.isdir(TARGET_ROI_PATH):
        print(f"Error: Target ROI directory {TARGET_ROI_PATH} not found. Exiting.")
        return

    # Define the full path for the output visualization (sub)folder
    # It will be inside TARGET_ROI_PATH/visualization/OUTPUT_VIS_SUBFOLDER
    vis_parent_dir = os.path.join(TARGET_ROI_PATH, "visualization") 
    final_output_dir = os.path.join(vis_parent_dir, OUTPUT_VIS_SUBFOLDER)

    try:
        os.makedirs(final_output_dir, exist_ok=True)
        print(f"Ensured output directory exists: {final_output_dir}")
    except OSError as e:
        print(f"Error: Could not create output directory {final_output_dir}. {e}")
        return

    visualize_single_roi(TARGET_ROI_PATH, final_output_dir)
    print("Visualization process completed.")

if __name__ == '__main__':
    main() 