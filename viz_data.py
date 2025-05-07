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
BASE_DATA_DIR = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_crawled" # Absolute path to data_crawled
OUTPUT_VIS_DIR = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/visualization_plots" # Absolute path for saving plots

# Keywords for dynamic folder searching or fixed names for direct path construction
NDVI_FOLDER_KEYWORD = "ndvi8days"
RVI_FOLDER_KEYWORD = "rvi_8days"
LST_FIXED_FOLDER_NAME = "lst"
ERA5_FIXED_FOLDER_NAME = "era5"
# MODIS_FIXED_FOLDER_NAME = "modis" # Commented out for MODIS data

ERA5_COMMON_PATTERN = "ERA5_temperature_skin_*.tif"
# MODIS_LST_DAY_PATTERN = "MODIS_LST_Day_*.tif" # Commented out for MODIS LST Day files

BANDS_CONFIG = {
    "RVI": {"folder_keyword": RVI_FOLDER_KEYWORD, "pattern": "*.tif", "band_index": 1},
    "NDVI": {"folder_keyword": NDVI_FOLDER_KEYWORD, "pattern": "*.tif", "band_index": 1},
    "LST": {"folder": LST_FIXED_FOLDER_NAME, "pattern": "lst16days_*.tif", "band_index": 1},
    # "MODIS_LST_Day": {"folder": MODIS_FIXED_FOLDER_NAME, "pattern": MODIS_LST_DAY_PATTERN, "band_index": 1}, # Commented out MODIS entry
    "ERA5_T2M": {"folder": ERA5_FIXED_FOLDER_NAME, "pattern": ERA5_COMMON_PATTERN, "band_index": 1},
    "ERA5_SKT": {"folder": ERA5_FIXED_FOLDER_NAME, "pattern": ERA5_COMMON_PATTERN, "band_index": 2},
}
BANDS_ORDER = ["RVI", "NDVI", "LST", "ERA5_T2M", "ERA5_SKT"] # Removed MODIS_LST_Day, back to 5 rows

NUM_RANDOM_ROIS = 3
NUM_IMAGES_PER_BAND = 10 # 10 columns

# --- Helper Functions ---

def extract_date_from_filename(filename):
    """Extracts date (YYYY-MM-DD) from common filename patterns."""
    basename = os.path.basename(filename)
    patterns = [
        r"(\d{4}-\d{2}-\d{2}).tif",
        r"_(\d{4}-\d{2}-\d{2}).tif",
        r"(\d{8}).tif",
        r"_(\d{8})_"
    ]
    for pattern in patterns:
        match = re.search(pattern, basename)
        if match:
            date_str = match.group(1)
            if len(date_str) == 8:
                return f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}"
            return date_str
    return os.path.splitext(basename)[0][:10]

def get_band_files(roi_path, band_name_key):
    """Gets all tif files for a specific band configuration in an ROI,
       dynamically finding the folder if a keyword is provided."""
    config = BANDS_CONFIG[band_name_key]
    actual_folder_path = None

    if "folder_keyword" in config:
        keyword = config["folder_keyword"].lower()
        try:
            found_dir = False
            for sub_dir_name in os.listdir(roi_path):
                sub_dir_full_path = os.path.join(roi_path, sub_dir_name)
                if os.path.isdir(sub_dir_full_path) and keyword in sub_dir_name.lower():
                    actual_folder_path = sub_dir_full_path
                    found_dir = True
                    break # Use the first directory found matching the keyword
            if not found_dir:
                # print(f"Info: No subfolder containing keyword '{keyword}' found in {roi_path} for {band_name_key}")
                return []
        except OSError as e:
            # print(f"OSError while searching for keyword '{keyword}' in {roi_path}: {e}")
            return []
    elif "folder" in config:
        fixed_folder_name = config["folder"]
        actual_folder_path = os.path.join(roi_path, fixed_folder_name)
        if not os.path.isdir(actual_folder_path):
            # print(f"Info: Fixed folder {actual_folder_path} not found for {band_name_key}")
            return []
    else:
        # print(f"Error: Band config for {band_name_key} lacks 'folder_keyword' or 'folder'.")
        return []

    search_pattern = os.path.join(actual_folder_path, config["pattern"])
    files = glob.glob(search_pattern)
    # print(f"Found {len(files)} files for {band_name_key} in {actual_folder_path} using pattern {config['pattern']}")
    return sorted(files)

def find_eligible_rois(base_dir):
    """Finds ROIs that have data files for all configured bands."""
    eligible_rois = []
    if not os.path.isdir(base_dir):
        print(f"Error: Base data directory {base_dir} not found.")
        return []

    for roi_name in os.listdir(base_dir):
        roi_path = os.path.join(base_dir, roi_name)
        if not os.path.isdir(roi_path):
            continue

        has_all_band_files = True
        for band_key in BANDS_ORDER: # Check based on BANDS_ORDER which implies specific band needs
            files = get_band_files(roi_path, band_key)
            if not files:
                has_all_band_files = False
                break
        
        if has_all_band_files:
            eligible_rois.append(roi_path)
            
    return eligible_rois

def visualize_roi(roi_path, vis_output_dir_roi):
    """Creates and saves a plot for a single ROI."""
    roi_name = os.path.basename(roi_path)
    print(f"  Visualizing ROI: {roi_name}")

    fig, axes = plt.subplots(len(BANDS_ORDER), NUM_IMAGES_PER_BAND, 
                             figsize=(NUM_IMAGES_PER_BAND * 3, len(BANDS_ORDER) * 3),
                             squeeze=False)
    
    plt.suptitle(f"Data Visualization for ROI: {roi_name}", fontsize=16, y=1.02)

    for i, band_name_key in enumerate(BANDS_ORDER):
        current_band_config = BANDS_CONFIG[band_name_key]
        band_files = get_band_files(roi_path, band_name_key) # Gets files matching the pattern for this band config
        
        if not band_files:
            axes[i, 0].set_ylabel(f"{band_name_key}\n(No data)", fontsize=10, rotation=0, labelpad=40, ha='right')
            for j in range(NUM_IMAGES_PER_BAND):
                axes[i, j].axis('off')
            continue

        axes[i, 0].set_ylabel(band_name_key, fontsize=10, rotation=0, labelpad=40, ha='right')

        if len(band_files) > NUM_IMAGES_PER_BAND:
            selected_files = random.sample(band_files, NUM_IMAGES_PER_BAND)
        else:
            selected_files = band_files
        
        selected_files.sort()

        for j in range(NUM_IMAGES_PER_BAND):
            ax = axes[i, j]
            ax.axis('off')

            if j < len(selected_files):
                file_path = selected_files[j]
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
                            vmin = np.nanpercentile(data, 2)
                            vmax = np.nanpercentile(data, 98)
                            if vmin == vmax:
                                vmin = np.nanmin(data)
                                vmax = np.nanmax(data)
                                if vmin == vmax:
                                     vmin = vmin - 0.1
                                     vmax = vmax + 0.1 if vmax != 0 else 0.1
                            
                            im = None # Initialize im
                            if vmin is None or vmax is None or vmin == vmax :
                                im = ax.imshow(data, cmap='viridis')
                            else:
                                im = ax.imshow(data, cmap='viridis', vmin=vmin, vmax=vmax)
                            
                            if im: # Add colorbar if image was plotted
                                plt.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)

                    file_date = extract_date_from_filename(file_path)
                    ax.set_title(file_date, fontsize=8)
                except IndexError:
                    ax.text(0.5, 0.5, f'Band {current_band_config["band_index"]}\nNot Found', ha='center', va='center', transform=ax.transAxes, fontsize=8)
                    print(f"    Band {current_band_config['band_index']} not found in {file_path}")
                except RasterioIOError:
                    ax.text(0.5, 0.5, 'Read Error', ha='center', va='center', transform=ax.transAxes)
                    print(f"    Error reading: {file_path}")
                except Exception as e:
                    ax.text(0.5, 0.5, 'Plot Error', ha='center', va='center', transform=ax.transAxes)
                    print(f"    Error plotting {file_path} for band {band_name_key}: {e}")
            else:
                 pass

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_png_path = os.path.join(vis_output_dir_roi, f"{roi_name}_visualization.png")
    try:
        plt.savefig(output_png_path, dpi=150, bbox_inches='tight')
        print(f"    Saved plot to: {output_png_path}")
    except Exception as e:
        print(f"    Error saving plot {output_png_path}: {e}")
    plt.close(fig)

def main():
    print("Starting visualization process...")
    print(f"Base data directory: {BASE_DATA_DIR}")
    print(f"Output directory for visualizations: {OUTPUT_VIS_DIR}")

    if os.path.exists(OUTPUT_VIS_DIR):
        print(f"Output directory {OUTPUT_VIS_DIR} already exists. Files might be overwritten.")
    else:
        try:
            os.makedirs(OUTPUT_VIS_DIR, exist_ok=True)
            print(f"Created output directory: {OUTPUT_VIS_DIR}")
        except OSError as e:
            print(f"Error: Could not create output directory {OUTPUT_VIS_DIR}. {e}")
            return

    eligible_rois = find_eligible_rois(BASE_DATA_DIR)

    if not eligible_rois:
        print("No ROIs found with all required data bands. Exiting.")
        return

    print(f"Found {len(eligible_rois)} eligible ROIs.")

    if len(eligible_rois) > NUM_RANDOM_ROIS:
        selected_rois = random.sample(eligible_rois, NUM_RANDOM_ROIS)
        print(f"Randomly selected {NUM_RANDOM_ROIS} ROIs for visualization.")
    else:
        selected_rois = eligible_rois
        print(f"Visualizing all {len(selected_rois)} eligible ROIs (fewer than or equal to {NUM_RANDOM_ROIS}).")

    for roi_path in selected_rois:
        roi_name = os.path.basename(roi_path)
        vis_output_dir_roi = os.path.join(OUTPUT_VIS_DIR, roi_name)
        try:
            os.makedirs(vis_output_dir_roi, exist_ok=True)
        except OSError as e:
            print(f"  Error creating directory {vis_output_dir_roi}. Skipping ROI {roi_name}. {e}")
            continue
        
        visualize_roi(roi_path, vis_output_dir_roi)

    print("Visualization process completed.")

if __name__ == '__main__':
    main()
