import os
import glob
import time
from datetime import datetime, timedelta
from pathlib import Path
import rasterio
from rasterio.warp import transform_bounds
from lst_retrieval import lst_retrive
import ee
# Note: Removed rvi_retrieval import as it's not used in the main flow requested
# from rvi_retrieval import get_sentinel1_collection, calculate_8day_composites_sar, sort_by_time, smooth_time_series
import shutil
# Initialize Earth Engine

X = 40
Y = 70
SKIPPING = True
TARGET_CRS = 'EPSG:4326' # Define target CRS globally
EXPORT_SCALE = 10 # Define export scale globally

try:
    ee.Initialize(project='ee-hadat')
except Exception:
    ee.Authenticate()
    ee.Initialize(project='ee-hadat')

def get_roi_coords_from_tif(tif_path):
    """Reads bounds from a TIF, converts to TARGET_CRS if needed."""
    with rasterio.open(tif_path) as dataset:
        bounds = dataset.bounds
        # Transform bounds to the globally defined TARGET_CRS if they are different
        if dataset.crs.to_string() != TARGET_CRS:
            print(f"Transforming bounds from {dataset.crs} to {TARGET_CRS}")
            bounds = transform_bounds(dataset.crs, TARGET_CRS, *bounds)
        # Coordinates format for ee.Geometry.Polygon: [[x, y], ...]
        coordinates = [
            [bounds[0], bounds[1]], # Lower Left
            [bounds[2], bounds[1]], # Lower Right
            [bounds[2], bounds[3]], # Upper Right
            [bounds[0], bounds[3]], # Upper Left
            [bounds[0], bounds[1]]  # Close the loop
        ]
        # Ensure coordinates are standard Python floats, not numpy types
        coordinates = [[float(x), float(y)] for x, y in coordinates]
        return coordinates

def get_ndvi_dates(ndvi_folder):
    """Gets sorted list of dates from NDVI filenames."""
    # Search for common NDVI naming patterns
    tif_files = sorted(glob.glob(os.path.join(ndvi_folder, 'ndvi*.tif')) + 
                       glob.glob(os.path.join(ndvi_folder, '*ndvi*.tif')))
    dates = []
    for tif in tif_files:
        base = os.path.basename(tif)
        try:
            # Attempt to extract date (assuming YYYY-MM-DD format at the end)
            date_str = base.split('_')[-1].replace('.tif', '')
            date = datetime.strptime(date_str, '%Y-%m-%d')
            dates.append(date)
        except ValueError: # Handle cases where parsing fails
            print(f"Could not parse date from filename: {base}")
            continue
        except Exception as e:
            print(f"Error processing file {base}: {e}")
            continue
    return sorted(list(set(dates))) # Use set to ensure uniqueness

def get_lst_dates(lst_folder):
    """Gets sorted list of dates from LST filenames."""
    # Search for common LST naming patterns (e.g., lst16days_YYYY-MM-DD.tif)
    tif_files = sorted(glob.glob(os.path.join(lst_folder, 'lst*.tif')))
    dates = []
    for tif in tif_files:
        base = os.path.basename(tif)
        try:
            # Attempt to extract date (assuming YYYY-MM-DD format at the end)
            date_str = base.split('_')[-1].replace('.tif', '')
            date = datetime.strptime(date_str, '%Y-%m-%d')
            dates.append(date)
        except ValueError: # Handle cases where parsing fails
            print(f"Could not parse date from LST filename: {base}")
            continue
        except Exception as e:
            print(f"Error processing LST file {base}: {e}")
            continue
    return sorted(list(set(dates))) # Use set to ensure uniqueness

def generate_missing_dates(date_list, max_interval=16, step_days=8):
    """Generate a list of dates to fill gaps larger than max_interval."""
    if not date_list or len(date_list) < 2:
        return []
    
    missing_dates = []
    date_list = sorted(list(set(date_list))) # Ensure sorted unique dates
    
    for i in range(len(date_list) - 1):
        current = date_list[i]
        next_date = date_list[i + 1]
        days_diff = (next_date - current).days
        
        # If the gap is large enough, add intermediate dates at step_days interval
        if days_diff >= max_interval:
            num_steps = days_diff // step_days
            for j in range(1, num_steps + 1):
                new_date = current + timedelta(days=j * step_days)
                if new_date < next_date: # Ensure we don't overshoot
                     # Check if new_date is too close to an existing date
                    is_close = False
                    for existing_date in date_list + missing_dates:
                        if abs((new_date - existing_date).days) <= step_days // 2:
                           is_close = True
                           break
                    if not is_close:
                           missing_dates.append(new_date)
    
    # Return unique sorted list of missing dates
    return sorted(list(set(missing_dates)))


def export_ee_image(image, bands, region, out_path, scale, crs=TARGET_CRS): # Default CRS set
    """Exports an EE image, handling potential errors and file types."""
    # Create proper geometry object with explicit projection (using the function's crs)
    try:
        region_geometry = ee.Geometry.Polygon(region, proj=crs, evenOdd=False) # Use proj parameter
    except Exception as e:
         print(f"Error creating geometry for {out_path}: {e}")
         print(f"Region data: {region}")
         print(f"CRS: {crs}")
         # Attempt with potential coordinate order swap for safety?
         # try:
         #      swapped_region = [[y, x] for x, y in region]
         #      region_geometry = ee.Geometry.Polygon(swapped_region, proj=crs, evenOdd=False)
         #      print("Used swapped coordinates for geometry")
         # except Exception as e2:
         #      print(f"Error creating geometry even with swapped coords: {e2}")
         #      return # Cannot proceed without valid geometry

         # Let's assume the original order is correct based on get_roi_coords_from_tif
         # and fail if it's wrong. Re-raise or return.
         return


    
    # Clip the image first
    image = image.clip(region_geometry)

    # Check if image has bands after clipping (can become empty)
    try:
        band_info = image.bandNames().getInfo()
        if not band_info or not any(b in band_info for b in bands):
             print(f"Warning: Image for {out_path} has no bands or missing required bands after clipping. Skipping export.")
             return
        # Select only available bands requested
        available_bands = [b for b in bands if b in band_info]
        if not available_bands:
             print(f"Warning: None of the requested bands {bands} are available in the clipped image for {out_path}. Skipping.")
             return
        image = image.select(available_bands)
    except ee.EEException as e:
         print(f"EE Error checking bands for {out_path}: {e}. Skipping export.")
         return
    except Exception as e:
         print(f"Unexpected error checking bands for {out_path}: {e}. Skipping export.")
         return


    
    # Then create download URL with proper parameters
    download_options = {
        'scale': scale,
        'region': region_geometry.getInfo()['coordinates'], # Pass coordinates explicitly
        'fileFormat': 'GeoTIFF',
        'maxPixels': 1e13,
        'format': 'GEO_TIFF', # Redundant? GeoTIFF specified above
        'crs': crs, # Explicitly set CRS
        'expires': 3600, # Increase expiry time?
    }

    try:
        url = image.getDownloadURL(download_options)
    except ee.EEException as e:
         print(f"EE Error getting download URL for {out_path}: {e}")
         # Print details that might help debugging EE issues
         print(f"Download Options: scale={scale}, crs={crs}, bands={available_bands}")
         # print(f"Region Coords: {region_geometry.getInfo()['coordinates']}") # Can be very large
         print(f"Image Info: {image.getInfo()}") # Can be large
         return
    except Exception as e:
         print(f"Unexpected error getting download URL for {out_path}: {e}")
         return

    
    import requests, zipfile, shutil # Keep imports local to function
    temp_dir = os.path.dirname(out_path) + '/temp_dl'
    os.makedirs(temp_dir, exist_ok=True)
    temp_file_path = os.path.join(temp_dir, 'download.tif')
    
    download_success = False
    try:
        print(f"Attempting download for {os.path.basename(out_path)}...")
        r = requests.get(url, stream=True, timeout=600, headers={'User-Agent': 'Mozilla/5.0'}) # Increased timeout
        
        # Check if the response is valid before saving
        if r.status_code != 200:
            # Try to get more info from the response if possible
            error_message = r.text[:500] if r.text else "(No response body)"
            raise Exception(f"Bad response: {r.status_code}, {error_message}")
            
        with open(temp_file_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024*1024): # Larger chunk size
                if chunk:
                    f.write(chunk)
        
        # Check if file is very small (potential error indicator)
        if os.path.getsize(temp_file_path) < 1024: # 1KB threshold
             print(f"Warning: Downloaded file {temp_file_path} is very small ({os.path.getsize(temp_file_path)} bytes). May be an error.")
             # Optionally read content if small
             with open(temp_file_path, 'r') as f_small:
                  print(f"Content of small file: {f_small.read(500)}")
             # Decide whether to proceed or raise error
             # raise Exception("Downloaded file too small")


        # Check file type by looking at its header
        with open(temp_file_path, 'rb') as f:
            header = f.read(4)
        
        # Handle different file types
        if zipfile.is_zipfile(temp_file_path):
            print("Downloaded a ZIP file unexpectedly, extracting...")
            with zipfile.ZipFile(temp_file_path, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            tif_files = [f for f in os.listdir(temp_dir) if f.endswith('.tif')]
            if tif_files:
                # Assuming the first TIF is the one we want
                shutil.move(os.path.join(temp_dir, tif_files[0]), out_path)
                print(f"Extracted and saved to {out_path}")
                download_success = True
            else:
                 print("Error: ZIP file did not contain a .tif file.")
        elif header[:2] in [b'II', b'MM']: # TIFF header magic numbers
            # It's a TIFF file
            shutil.move(temp_file_path, out_path)
            print(f"Saved direct TIFF to {out_path}")
            download_success = True
        else:
            # Unknown file type
            print(f"Downloaded file is not recognized as TIFF or ZIP: {temp_file_path}")
            with open(temp_file_path, 'rb') as f:
                print(f"First 100 bytes: {f.read(100)}")
                
    except requests.exceptions.RequestException as e:
         print(f"Download Error (requests): {out_path}: {e}")
    except Exception as e:
        print(f"Error processing download for {out_path}: {e}")
    finally:
        # Clean up temp dir
        if os.path.exists(temp_dir):
             shutil.rmtree(temp_dir, ignore_errors=True)
        # Optional: Verify image if download seemed successful
        # if download_success:
        #      verify_image(out_path) # Be cautious, this reads the file again


def verify_image(img_path):
    """Verify the image has been properly saved and print its CRS and shape"""
    try:
        with rasterio.open(img_path) as src:
            crs = src.crs
            shape = src.shape
            bounds = src.bounds
            if not crs:
                #  print(f"Image verification - {os.path.basename(img_path)}:")
                 print(f"  Warning: CRS is missing!")
                #  print(f"  Shape: {shape}")
                #  print(f"  Bounds: {bounds}")
                 return False # Indicate missing CRS
            # print(f"Image verification - {os.path.basename(img_path)}:")
            print(f"  CRS: {crs}")
            # print(f"  Shape: {shape}")
            # print(f"  Bounds: {bounds}")
            return True # Indicate success
    except rasterio.errors.RasterioIOError as e:
         print(f"Rasterio Error verifying image {img_path}: {e}")
         return False
    except Exception as e:
        print(f"Unexpected Error verifying image {img_path}: {e}")
        return False

def reproject_if_needed(src_path, target_crs=TARGET_CRS): # Default target_crs
    """Reproject an image to the target CRS if needed"""
    try:
        with rasterio.open(src_path) as src:
            # Check if CRS is missing or different
            if not src.crs or src.crs != target_crs:
                action = "Assigning" if not src.crs else "Reprojecting"
                source_crs_str = src.crs if src.crs else "Missing CRS"
                print(f"{action} {os.path.basename(src_path)} from {source_crs_str} to {target_crs}")
                
                # Create output path
                dst_path = src_path.replace('.tif', '_reprojected.tif')
                
                # Get transformation parameters
                # Handle missing source CRS - default to assuming TARGET_CRS? Or fail?
                # For now, let's assume if missing, we just assign. Reprojection needs source.
                if not src.crs:
                     print("Warning: Source CRS missing. Assigning target CRS without reprojection.")
                     # Create a new profile with the target CRS
                     dst_kwargs = src.meta.copy()
                     dst_kwargs['crs'] = target_crs
                     # Read data and write to new file with assigned CRS
                     data = src.read()
                     with rasterio.open(dst_path, 'w', **dst_kwargs) as dst:
                          dst.write(data)
                else:
                     # Proceed with reprojection
                     transform, width, height = rasterio.warp.calculate_default_transform(
                         src.crs, target_crs, src.width, src.height, *src.bounds
                     )
                     
                     # Set metadata for output
                     dst_kwargs = src.meta.copy()
                     dst_kwargs.update({
                         'crs': target_crs,
                         'transform': transform,
                         'width': width,
                         'height': height,
                         'nodata': src.nodata # Preserve nodata value if exists
                     })
                     
                     # Reproject and save
                     with rasterio.open(dst_path, 'w', **dst_kwargs) as dst:
                         for i in range(1, src.count + 1):
                             rasterio.warp.reproject(
                                 source=rasterio.band(src, i),
                                 destination=rasterio.band(dst, i),
                                 src_transform=src.transform,
                                 src_crs=src.crs,
                                 src_nodata=src.nodata, # Pass nodata for source
                                 dst_transform=transform,
                                 dst_crs=target_crs,
                                 dst_nodata=dst.nodata, # Pass nodata for destination
                                 resampling=rasterio.warp.Resampling.nearest # Use nearest for categorical/thematic, bilinear for continuous
                                                                             # Let's default to nearest for safety, maybe bilinear for LST/Temp?
                             )
                
                # Replace original with reprojected/assigned version
                shutil.move(dst_path, src_path)
                print(f"Successfully updated CRS for {os.path.basename(src_path)}")
                return True # Indicated a change was made
            return False  # No reprojection needed
    except rasterio.errors.RasterioIOError as e:
        print(f"Rasterio Error processing reprojection for {src_path}: {e}")
        return False # Indicate error/no change
    except Exception as e:
        print(f"Unexpected Error reprojecting {src_path}: {e}")
        return False # Indicate error/no change

def get_era5_for_date(target_date, roi_geom, region, out_folder, window=8):
    """Fetches and exports the closest ERA5 Land image for a target date."""
    start = ee.Date((target_date - timedelta(days=window)).strftime('%Y-%m-%d'))
    end = ee.Date((target_date + timedelta(days=window)).strftime('%Y-%m-%d'))
    
    try:
        era5 = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR') \
            .filterDate(start, end) \
            .filterBounds(roi_geom) # Filter by geometry first

        # Check if collection is empty
        collection_size = era5.size().getInfo()
        if collection_size == 0:
            print(f"No ERA5 images found for {target_date.strftime('%Y-%m-%d')} window.")
            return

        # Convert target_date to milliseconds since epoch for EE
        target_millis = ee.Date(target_date.strftime('%Y-%m-%d')).millis()
        
        # Define sort function that works with Earth Engine serialization
        def add_diff(img):
            diff = ee.Number(img.get('system:time_start')).subtract(target_millis).abs()
            return img.set('date_diff', diff)
        
        # Apply the function to add the date_diff property to each image
        era5_with_diff = era5.map(add_diff)
        
        # Sort by the date_diff property
        era5_sorted = era5_with_diff.sort('date_diff')
        best_img = ee.Image(era5_sorted.first())
        
        # Export the best image
        date_str = target_date.strftime('%Y-%m-%d')
        out_path = os.path.join(out_folder, f'ERA5_temperature_skin_{date_str}.tif')
        
        if SKIPPING and os.path.exists(out_path):
            print(f"ERA5 file for {date_str} already exists. Skipping download.")
            verify_image(out_path) # Optional: Verify existing file
            return
            
        print(f"Exporting ERA5 for {date_str}...")
        # Use globally defined scale and target CRS
        export_ee_image(best_img, ['temperature_2m', 'skin_temperature'], region, out_path, scale=EXPORT_SCALE, crs=TARGET_CRS)
        time.sleep(1) # Pause after export attempt

    except ee.EEException as e:
         print(f"EE Error fetching ERA5 for {target_date.strftime('%Y-%m-%d')}: {e}")
    except Exception as e:
         print(f"Unexpected error fetching ERA5 for {target_date.strftime('%Y-%m-%d')}: {e}")


def mask_modis_lst(img, day=True):
    """Masks MODIS LST based on QC flags."""
    if day:
        lst_band = 'LST_Day_1km'
        qc_band = 'QC_Day'
    else:
        lst_band = 'LST_Night_1km'
        qc_band = 'QC_Night'
        
    if lst_band not in img.bandNames().getInfo():
        # print(f"Warning: Band {lst_band} not found in MODIS image.")
        return None # Or return an empty image? ee.Image().select([])

    lst = img.select(lst_band)
    
    if qc_band not in img.bandNames().getInfo():
        print(f"Warning: Band {qc_band} not found in MODIS image. Returning LST unmasked.")
        return lst # Return LST without masking if QC is missing

    qc = img.select(qc_band)
    # Bitmask check for Mandatory QA flags (Bits 0-1) == 0 ('LST produced, good quality')
    mask = qc.bitwiseAnd(3).eq(0)
    return lst.updateMask(mask)


def get_modis_for_date(target_date, roi_geom, region, out_folder, window=8):
    """Fetches and exports the closest MODIS LST images for a target date."""
    start = ee.Date((target_date - timedelta(days=window)).strftime('%Y-%m-%d'))
    end = ee.Date((target_date + timedelta(days=window)).strftime('%Y-%m-%d'))
    
    date_str = target_date.strftime('%Y-%m-%d')
    out_path_day = os.path.join(out_folder, f'MODIS_LST_Day_{date_str}.tif')
    out_path_night = os.path.join(out_folder, f'MODIS_LST_Night_{date_str}.tif')
    
    if SKIPPING and os.path.exists(out_path_day) and os.path.exists(out_path_night):
        print(f"MODIS files for {date_str} already exist. Skipping download.")
        # Optional: Verify existing files
        # verify_image(out_path_day)
        # verify_image(out_path_night)
        return
        
    try:
        modis = ee.ImageCollection('MODIS/061/MOD11A1') \
            .filterDate(start, end) \
            .filterBounds(roi_geom) # Filter by geometry first

        # Check if collection is empty
        collection_size = modis.size().getInfo()
        if collection_size == 0:
            print(f"No MODIS images found for {date_str} window.")
            return
            
        modis_list = modis.toList(collection_size)
        best_img_info = None
        min_diff = float('inf') # Use float infinity

        # Iterate on the server side using map and aggregate
        def find_closest(image, target_millis_):
             img_millis = ee.Number(image.get('system:time_start'))
             diff = img_millis.subtract(target_millis_).abs()
             return image.set('date_diff', diff)

        target_millis = ee.Date(date_str).millis()
        modis_with_diff = modis.map(lambda img: find_closest(img, target_millis))
        best_img = ee.Image(modis_with_diff.sort('date_diff').first())


        # Check if a best image was found
        if not best_img.bandNames().size().getInfo() > 0: # Check if image has bands
             print(f"Could not find a valid best MODIS image for {date_str}")
             return

        print(f"Best MODIS image found for {date_str}: {best_img.id().getInfo()}")
        
        lst_day = mask_modis_lst(best_img, day=True)
        lst_night = mask_modis_lst(best_img, day=False)
        
        # Export Day LST if available
        if lst_day:
            day_bands = lst_day.bandNames().getInfo()
            if day_bands: # Ensure lst_day is not empty
                 print(f"Exporting MODIS Day LST for {date_str}...")
                 export_ee_image(lst_day, day_bands, region, out_path_day, scale=EXPORT_SCALE, crs=TARGET_CRS)
                 time.sleep(1)
            else:
                 print(f"MODIS Day LST masked out or unavailable for {date_str}.")
        else:
             print(f"MODIS Day LST band not found or masking failed for {date_str}.")

        # Export Night LST if available
        if lst_night:
            night_bands = lst_night.bandNames().getInfo()
            if night_bands: # Ensure lst_night is not empty
                 print(f"Exporting MODIS Night LST for {date_str}...")
                 export_ee_image(lst_night, night_bands, region, out_path_night, scale=EXPORT_SCALE, crs=TARGET_CRS)
                 time.sleep(1)
            else:
                 print(f"MODIS Night LST masked out or unavailable for {date_str}.")
        else:
             print(f"MODIS Night LST band not found or masking failed for {date_str}.")

    except ee.EEException as e:
         print(f"EE Error fetching MODIS for {date_str}: {e}")
    except Exception as e:
         print(f"Unexpected error fetching MODIS for {date_str}: {e}")


def main():
    base_dir = '/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_crawled'
    export_base = '/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_crawled'
    os.makedirs(export_base, exist_ok=True)
    date_start = '2022-12-15'
    date_end = '2025-01-01'


    # Get list of directories to process
    try:
        all_roi_folders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    except OSError as e:
         print(f"Error listing base directory {base_dir}: {e}")
         return

    # Process selected range [X:Y]
    roi_folders_to_process = all_roi_folders[X:Y]
    print(f"Processing {len(roi_folders_to_process)} ROIs from index {X} to {Y-1}")

    for roi_folder in roi_folders_to_process:
        print(f"\n===== Processing ROI: {roi_folder} =====")
        roi_path = os.path.join(base_dir, roi_folder)
        
        # Find NDVI directory (look for 'ndvi' in name)
        ndvi_dir = None
        potential_ndvi_dirs = [d for d in os.listdir(roi_path) if os.path.isdir(os.path.join(roi_path, d)) and 'ndvi8days' in d.lower()] # Looser check for 'ndvi'
        if not potential_ndvi_dirs:
            print(f"Warning: No NDVI directory found in {roi_path}. Skipping ROI.")
            continue
        ndvi_dir = os.path.join(roi_path, potential_ndvi_dirs[0]) # Use first found NDVI dir
        
        # Find any TIF file in the NDVI dir to get coordinates
        ndvi_tif_files = list(Path(ndvi_dir).glob('*.tif'))
        if not ndvi_tif_files:
            print(f"Warning: No .tif files found in {ndvi_dir}. Cannot determine ROI. Skipping ROI.")
            continue
        
        # Get coordinates from NDVI file and define geometry using the target CRS
        try:
            # Use the first TIF found in the NDVI directory
            coords = get_roi_coords_from_tif(str(ndvi_tif_files[0])) 
            # Define EE Geometry using the coordinates and specifying they are in TARGET_CRS
            roi_geom = ee.Geometry.Polygon(coords, proj=TARGET_CRS, evenOdd=False)
            # Region for export can be the raw coordinates list
            region = coords
            print(f"ROI Geometry defined using {TARGET_CRS} coordinates from {os.path.basename(ndvi_tif_files[0])}")
        except Exception as e:
            print(f"Error getting coordinates or defining geometry for {roi_folder}: {e}. Skipping ROI.")
            continue
        
        # out_folder_lst = os.path.join(export_base, roi_folder, 'lst')
        # os.makedirs(out_folder_lst, exist_ok=True)
        
        lst_retrive(date_start, date_end,roi_geom, roi_folder, base_dir, SKIPPING)
        
        # --- Find LST directory for dates ---
        lst_dir = None
        potential_lst_dirs = [d for d in os.listdir(roi_path) if os.path.isdir(os.path.join(roi_path, d)) and 'lst' in d.lower()] # Look for 'lst' dir
        if not potential_lst_dirs:
            print(f"Warning: No LST directory found in {roi_path} to determine dates. Skipping ROI.")
            continue
        lst_dir = os.path.join(roi_path, potential_lst_dirs[0]) # Use first found LST dir
        print(f"Using NDVI directory for ROI: {ndvi_dir}")
        print(f"Using LST directory for Dates: {lst_dir}")
            
        # Get original LST dates
        lst_dates = get_lst_dates(str(lst_dir)) # Get dates from LST directory
        if not lst_dates:
             print(f"Warning: No valid dates found from LST files in {lst_dir}. Skipping data download.")
             continue
        print(f"Found {len(lst_dates)} unique LST dates for {roi_folder}")
        
        # Generate missing dates between LST intervals
        missing_dates = generate_missing_dates(lst_dates) # Generate gaps based on LST dates
        print(f"Identified {len(missing_dates)} missing dates to fill gaps")
        
        # Combine all dates to process
        all_dates = sorted(list(set(lst_dates + missing_dates)))
        print(f"Total dates to process (original LST + gaps): {len(all_dates)}")
        
        out_folder_era5 = os.path.join(export_base, roi_folder, 'era5')
        # out_folder_modis = os.path.join(export_base, roi_folder, 'modis')
        os.makedirs(out_folder_era5, exist_ok=True)
        # os.makedirs(out_folder_modis, exist_ok=True)
        
        for date in all_dates:
            date_str = date.strftime('%Y-%m-%d')
            if date in lst_dates: # Check against original LST dates
                print(f'\n--- Processing Original LST date: {roi_folder}, Date: {date_str} ---')
            else:
                print(f'\n--- Processing Missing date (gap-fill): {roi_folder}, Date: {date_str} ---')
            
            # print(f'Fetching ERA5 for date: {date_str}')
            get_era5_for_date(date, roi_geom, region, out_folder_era5)
            
            # print(f'Fetching MODIS for date: {date_str}')
            # get_modis_for_date(date, roi_geom, region, out_folder_modis)
        
        # --- Post-processing: Check/Fix CRS ---
        print(f"\nVerifying/Fixing coordinate systems to {TARGET_CRS} for {roi_folder}...")
        
        # Check ERA5 files
        era5_files = glob.glob(os.path.join(out_folder_era5, '*.tif'))
        print(f"Checking {len(era5_files)} ERA5 files...")
        for era5_file in era5_files:
            reproject_if_needed(era5_file, TARGET_CRS)
            
        # Check MODIS files
        # modis_files = glob.glob(os.path.join(out_folder_modis, '*.tif'))
        # print(f"Checking {len(modis_files)} MODIS files...")
        # for modis_file in modis_files:
        #     reproject_if_needed(modis_file, TARGET_CRS)
        
        print(f"===== Finished processing ROI: {roi_folder} =====")
        print("--------------------------------")

if __name__ == '__main__':
    main()
