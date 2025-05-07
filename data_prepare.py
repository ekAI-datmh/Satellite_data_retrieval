import os
import shutil
import glob
from datetime import datetime, timedelta, date
import rasterio
import numpy as np
import re

# --- Configuration ---
SOURCE_BASE_DIR = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_crawled"
TARGET_BASE_DIR = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_lst_16days"

LST_SRC_FOLDER_NAME = "lst"  # Source LST folder name in data_crawled/ROI/
ERA5_SRC_FOLDER_NAME = "era5" # Source ERA5 folder name in data_crawled/ROI/

# Target folder names within data_lst_16days/ROI/
LST_TRG_FOLDER_NAME = "lst"
ERA5_TRG_FOLDER_NAME = "era5"

START_DATE_SERIES = date(2023, 1, 1)
END_DATE_SERIES = date(2025, 1, 1)
TIME_INTERVAL_DAYS = 16
SEARCH_WINDOW_DAYS = 8  # +/- days around the target date

# --- Helper Functions ---

def parse_date_from_filename(filename_str, expected_prefix):
    """
    Parses YYYY-MM-DD from filenames like 'prefix_YYYY-MM-DD.tif'.
    Returns a datetime.date object or None.
    """
    # Strip prefix and .tif extension
    # Example: lst16days_2023-01-15.tif -> 2023-01-15
    # Example: ERA5_temperature_skin_2023-01-15.tif -> 2023-01-15 (needs more specific regex)
    
    if expected_prefix == "lst16days_":
        # Pattern for "lst16days_YYYY-MM-DD.tif"
        match = re.search(r"lst16days_(\d{4}-\d{2}-\d{2})\.tif$", filename_str)
        if match:
            date_part = match.group(1)
            try:
                return datetime.strptime(date_part, "%Y-%m-%d").date()
            except ValueError:
                return None
    elif expected_prefix == "ERA5_temperature_skin_":
         # Pattern for "ERA5_temperature_skin_YYYY-MM-DD.tif"
        match = re.search(r"ERA5_temperature_skin_(\d{4}-\d{2}-\d{2})\.tif$", filename_str)
        if match:
            date_part = match.group(1)
            try:
                return datetime.strptime(date_part, "%Y-%m-%d").date()
            except ValueError:
                return None
    return None


def calculate_nan_percentage(image_path):
    """
    Calculates the percentage of NaN or NoData pixels in a raster image.
    Returns percentage (0-100) or 101.0 if an error occurs.
    """
    try:
        with rasterio.open(image_path) as src:
            # Read the first band as a masked array.
            # The mask is True for invalid data (NoData or NaN).
            data_masked_array = src.read(1, masked=True)
            
            # Count masked (invalid/NoData/NaN) pixels
            nodata_pixels = np.sum(data_masked_array.mask)
            total_pixels = data_masked_array.size
            
            if total_pixels == 0:
                # print(f"Warning: Image {image_path} has zero total pixels.")
                return 101.0  # Indicate error or empty raster
            
            percentage = (nodata_pixels / total_pixels) * 100.0
            return percentage
    except Exception as e:
        print(f"Error calculating NaN for {image_path}: {e}")
        return 101.0  # Indicate error

def main_process():
    """
    Main processing function.
    """
    print(f"Starting data preparation process...")
    print(f"Source base directory: {SOURCE_BASE_DIR}")
    print(f"Target base directory: {TARGET_BASE_DIR}")

    if not os.path.isdir(SOURCE_BASE_DIR):
        print(f"Error: Source directory {SOURCE_BASE_DIR} not found.")
        return
    if not os.path.isdir(TARGET_BASE_DIR):
        print(f"Error: Target directory {TARGET_BASE_DIR} not found.")
        return

    roi_names_in_target = list(os.listdir(TARGET_BASE_DIR))[:1]

    for roi_name in roi_names_in_target:
        target_roi_path = os.path.join(TARGET_BASE_DIR, roi_name)
        
        # Skip files like .json, process only directories
        if not os.path.isdir(target_roi_path) or roi_name.lower().endswith('.json'):
            # print(f"Skipping non-directory or JSON file: {roi_name}")
            continue

        print(f"\nProcessing ROI: {roi_name}")

        source_roi_path = os.path.join(SOURCE_BASE_DIR, roi_name)
        if not os.path.isdir(source_roi_path):
            print(f"  Warning: Source ROI folder {source_roi_path} not found. Skipping this ROI.")
            continue

        # 1. Clean and create target LST and ERA5 folders
        target_lst_folder = os.path.join(target_roi_path, LST_TRG_FOLDER_NAME)
        target_era5_folder = os.path.join(target_roi_path, ERA5_TRG_FOLDER_NAME)

        for folder_to_recreate in [target_lst_folder, target_era5_folder]:
            if os.path.exists(folder_to_recreate):
                print(f"  Removing existing folder: {folder_to_recreate}")
                shutil.rmtree(folder_to_recreate)
            os.makedirs(folder_to_recreate)
            print(f"  Created empty folder: {folder_to_recreate}")
        
        # Pre-fetch and parse dates for all available LST files in the source ROI
        source_lst_dir = os.path.join(source_roi_path, LST_SRC_FOLDER_NAME)
        if not os.path.isdir(source_lst_dir):
            print(f"  Warning: Source LST directory {source_lst_dir} not found. Skipping LST processing for this ROI.")
            available_lst_files_info = []
        else:
            available_lst_files_info = []
            for f_name in os.listdir(source_lst_dir):
                if f_name.endswith(".tif"): # General check, specific prefix check done by parser
                    original_lst_date = parse_date_from_filename(f_name, "lst16days_")
                    if original_lst_date:
                        available_lst_files_info.append({
                            "path": os.path.join(source_lst_dir, f_name),
                            "original_date": original_lst_date
                        })
            if not available_lst_files_info:
                 print(f"  Warning: No valid LST files found in {source_lst_dir}.")


        # 2. Generate target time series and process each step
        current_target_date = START_DATE_SERIES
        while current_target_date <= END_DATE_SERIES:
            target_date_str = current_target_date.strftime("%Y-%m-%d")
            print(f"  Processing target date: {target_date_str}")

            window_start = current_target_date - timedelta(days=SEARCH_WINDOW_DAYS)
            window_end = current_target_date + timedelta(days=SEARCH_WINDOW_DAYS)
            
            candidate_lst_images = []
            if not available_lst_files_info: # Skip if no LST source files for ROI
                 print(f"    No source LST files for ROI {roi_name} to select from.")
                 current_target_date += timedelta(days=TIME_INTERVAL_DAYS)
                 continue


            for lst_info in available_lst_files_info:
                if window_start <= lst_info["original_date"] <= window_end:
                    day_diff = abs((lst_info["original_date"] - current_target_date).days)
                    
                    # print(f"      Considering LST: {lst_info['path']}, Date: {lst_info['original_date']}, DayDiff: {day_diff}")
                    nan_perc = calculate_nan_percentage(lst_info["path"])
                    
                    if nan_perc > 100.0: # Error in calculation or unreadable file
                        # print(f"        Skipping LST due to NaN calculation error: {lst_info['path']}")
                        continue
                        
                    score = day_diff * nan_perc # If nan_perc is 0, score is 0
                                                # If day_diff is 0 and nan_perc is low, score is low
                    candidate_lst_images.append({
                        "path": lst_info["path"],
                        "original_date": lst_info["original_date"],
                        "day_diff": day_diff,
                        "per_nan": nan_perc,
                        "score": score
                    })
            
            if not candidate_lst_images:
                print(f"    Warning: No suitable LST images found in window for target date {target_date_str}.")
                current_target_date += timedelta(days=TIME_INTERVAL_DAYS)
                continue

            # Select the best LST candidate (minimum score)
            # Tie-breaking: prefer smaller day_diff, then smaller per_nan
            best_lst_candidate = min(candidate_lst_images, key=lambda x: (x["score"], x["day_diff"], x["per_nan"]))
            print(f"    Selected LST: {os.path.basename(best_lst_candidate['path'])} (OrigDate: {best_lst_candidate['original_date']}, Score: {best_lst_candidate['score']:.2f})")

            # Copy and rename selected LST image
            target_lst_filename = f"lst_{target_date_str}.tif"
            destination_lst_path = os.path.join(target_lst_folder, target_lst_filename)
            try:
                shutil.copy2(best_lst_candidate["path"], destination_lst_path)
                # print(f"      Copied LST to: {destination_lst_path}")
            except Exception as e:
                print(f"      Error copying LST file {best_lst_candidate['path']}: {e}")
                current_target_date += timedelta(days=TIME_INTERVAL_DAYS)
                continue # Skip to next target date if LST copy fails

            # Find, copy, and rename corresponding ERA5 image
            lst_original_date_for_era5_search = best_lst_candidate["original_date"]
            source_era5_dir = os.path.join(source_roi_path, ERA5_SRC_FOLDER_NAME)
            found_era5_match = False

            if not os.path.isdir(source_era5_dir):
                print(f"    Warning: Source ERA5 directory {source_era5_dir} not found.")
            else:
                for era5_filename in os.listdir(source_era5_dir):
                    if era5_filename.endswith(".tif"): # General check
                        era5_file_date = parse_date_from_filename(era5_filename, "ERA5_temperature_skin_")
                        if era5_file_date and era5_file_date == lst_original_date_for_era5_search:
                            source_era5_path = os.path.join(source_era5_dir, era5_filename)
                            target_era5_filename = f"era5_{target_date_str}.tif"
                            destination_era5_path = os.path.join(target_era5_folder, target_era5_filename)
                            try:
                                shutil.copy2(source_era5_path, destination_era5_path)
                                # print(f"      Copied ERA5 to: {destination_era5_path}")
                                found_era5_match = True
                                break 
                            except Exception as e:
                                print(f"      Error copying ERA5 file {source_era5_path}: {e}")
                                break # Stop trying for this ERA5 if copy fails

            if not found_era5_match:
                print(f"    Warning: No corresponding ERA5 image found for LST's original date {lst_original_date_for_era5_search} (target date {target_date_str}).")

            current_target_date += timedelta(days=TIME_INTERVAL_DAYS)
            
    print("\nData preparation process completed.")

if __name__ == "__main__":
    main_process()
