import os
import glob
import time
from datetime import datetime, timedelta
from pathlib import Path
import rasterio
from rasterio.warp import transform_bounds
import ee
from rvi_retrieval import get_sentinel1_collection, calculate_8day_composites_sar, sort_by_time, smooth_time_series
import shutil
# Initialize Earth Engine

X = 27
Y = 53
try:
    ee.Initialize(project='ee-hadat')
except Exception:
    ee.Authenticate()
    ee.Initialize(project='ee-hadat')

def get_roi_coords_from_tif(tif_path):
    with rasterio.open(tif_path) as dataset:
        bounds = dataset.bounds
        if dataset.crs.to_string() != 'EPSG:3857':
            bounds = transform_bounds(dataset.crs, 'EPSG:3857', *bounds)
        coordinates = [
            [bounds[2], bounds[1]],
            [bounds[2], bounds[3]],
            [bounds[0], bounds[3]],
            [bounds[0], bounds[1]],
            [bounds[2], bounds[1]]
        ]
        return coordinates

def get_lst_dates(lst_folder):
    tif_files = sorted(glob.glob(os.path.join(lst_folder, '*.tif')))
    dates = []
    for tif in tif_files:
        base = os.path.basename(tif)
        try:
            date_str = base.split('_')[-1].replace('.tif', '')
            date = datetime.strptime(date_str, '%Y-%m-%d')
            dates.append(date)
        except Exception:
            continue
    return sorted(dates)

def generate_missing_dates(date_list, max_interval=16):
    """Generate a list of dates that are missing between the input dates.
    Only adds dates if the gap is larger than 8 days (or other threshold)."""
    if not date_list or len(date_list) < 2:
        return []
    
    missing_dates = []
    for i in range(len(date_list) - 1):
        current = date_list[i]
        next_date = date_list[i + 1]
        days_diff = (next_date - current).days
        
        # If the gap is large (e.g., > 8 days), add intermediate dates
        if days_diff >= 16:
            # Calculate how many intermediate points we want
            num_points = days_diff // 8  # This will divide the interval roughly into 8-day chunks
            
            if num_points > 0:
                interval = days_diff / (num_points)  # +1 because we already have the endpoints
                for j in range(1, num_points):
                    new_date = current + timedelta(days=int(j * interval))
                    missing_dates.append(new_date)
    
    return missing_dates

def export_ee_image(image, bands, region, out_path, scale, crs='EPSG:3857'):
    # Create proper geometry object with explicit projection
    region_geometry = ee.Geometry.Polygon([region], crs)
    
    # Clip the image first and make sure it has data
    image = image.clip(region_geometry)
    
    # Then create download URL with proper parameters
    url = image.select(bands).getDownloadURL({
        'scale': scale,
        'region': region_geometry,
        'fileFormat': 'GeoTIFF',  # Changed to GeoTIFF instead of ZIP
        'maxPixels': 1e13,
        'format': 'GEO_TIFF',
        'crs': crs,  # Explicitly set CRS
        'expires': 3600,
    })
    
    import requests, zipfile, shutil
    temp_dir = os.path.dirname(out_path) + '/temp_dl'
    os.makedirs(temp_dir, exist_ok=True)
    temp_file_path = os.path.join(temp_dir, 'download.tif')  # Changed from .zip to .tif
    try:
        r = requests.get(url, stream=True, timeout=300, headers={'User-Agent': 'Mozilla/5.0'})
        # Check if the response is valid before saving
        if r.status_code != 200:
            raise Exception(f"Bad response: {r.status_code}, {r.text}")
            
        with open(temp_file_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        # Check file type by looking at its header
        with open(temp_file_path, 'rb') as f:
            header = f.read(4)
        
        # Handle different file types
        if zipfile.is_zipfile(temp_file_path):
            # It's a ZIP file
            with zipfile.ZipFile(temp_file_path, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            tif_files = [f for f in os.listdir(temp_dir) if f.endswith('.tif')]
            if tif_files:
                shutil.move(os.path.join(temp_dir, tif_files[0]), out_path)
                print(f"Extracted and saved to {out_path}")
                # verify_image(out_path)
        elif header[:2] in [b'II', b'MM']:  # TIFF header magic numbers
            # It's a TIFF file
            shutil.move(temp_file_path, out_path)
            print(f"Saved direct TIFF to {out_path}")
            # verify_image(out_path)
        else:
            # Unknown file type
            print(f"Downloaded file is not recognized: {temp_file_path}")
            with open(temp_file_path, 'rb') as f:
                print(f"First 100 bytes: {f.read(100)}")
    except Exception as e:
        print(f"Error downloading/exporting {out_path}: {e}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def verify_image(img_path):
    """Verify the image has been properly saved and print its CRS and shape"""
    try:
        with rasterio.open(img_path) as src:
            crs = src.crs
            shape = src.shape
            bounds = src.bounds
            print(f"Image verification - {os.path.basename(img_path)}:")
            print(f"  CRS: {crs}")
            print(f"  Shape: {shape}")
            print(f"  Bounds: {bounds}")
    except Exception as e:
        print(f"Error verifying image {img_path}: {e}")

def reproject_if_needed(src_path, target_crs='EPSG:3857'):
    """Reproject an image to the target CRS if needed"""
    try:
        with rasterio.open(src_path) as src:
            if src.crs != target_crs:
                print(f"Reprojecting {os.path.basename(src_path)} from {src.crs} to {target_crs}")
                
                # Create output path
                dst_path = src_path.replace('.tif', '_reprojected.tif')
                
                # Get transformation parameters
                transform, width, height = rasterio.warp.calculate_default_transform(
                    src.crs, target_crs, src.width, src.height, *src.bounds
                )
                
                # Set metadata for output
                dst_kwargs = src.meta.copy()
                dst_kwargs.update({
                    'crs': target_crs,
                    'transform': transform,
                    'width': width,
                    'height': height
                })
                
                # Reproject and save
                with rasterio.open(dst_path, 'w', **dst_kwargs) as dst:
                    for i in range(1, src.count + 1):
                        rasterio.warp.reproject(
                            source=rasterio.band(src, i),
                            destination=rasterio.band(dst, i),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=transform,
                            dst_crs=target_crs,
                            resampling=rasterio.warp.Resampling.bilinear
                        )
                
                # Replace original with reprojected version
                shutil.move(dst_path, src_path)
                print(f"Successfully reprojected {os.path.basename(src_path)}")
                return True
            return False  # No reprojection needed
    except Exception as e:
        print(f"Error reprojecting {src_path}: {e}")
        return False

def get_era5_for_date(target_date, roi_geom, region, out_folder, window=8):
    start = ee.Date((target_date - timedelta(days=window)).strftime('%Y-%m-%d'))
    end = ee.Date((target_date + timedelta(days=window)).strftime('%Y-%m-%d'))
    era5 = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR') \
        .filterDate(start, end) \
        .filterBounds(roi_geom)
    
    # Convert target_date to milliseconds since epoch for EE
    target_millis = ee.Date(target_date.strftime('%Y-%m-%d')).millis()
    
    # Define sort function that works with Earth Engine serialization
    def add_diff(img):
        diff = ee.Number(ee.Date(img.get('system:time_start')).millis()).subtract(target_millis).abs()
        return img.set('date_diff', diff)
    
    # Apply the function to add the date_diff property to each image
    era5_with_diff = era5.map(add_diff)
    
    # Sort by the date_diff property
    era5_sorted = era5_with_diff.sort('date_diff')
    best_img = ee.Image(era5_sorted.first())
    
    # Export the best image
    date_str = target_date.strftime('%Y-%m-%d')
    out_path = os.path.join(out_folder, f'ERA5_temperature_skin_{date_str}.tif')
    
    # Skip if file already exists
    if os.path.exists(out_path):
        print(f"ERA5 file for {date_str} already exists. Skipping.")
        # verify_image(out_path)
        return
        
    # Use consistent CRS (EPSG:3857) for all exports
    export_ee_image(best_img, ['temperature_2m', 'skin_temperature'], region, out_path, 11132, crs='EPSG:3857')
    time.sleep(2)

def mask_modis_lst(img, day=True):
    if day:
        qc = img.select('QC_Day')
        mask = qc.bitwiseAnd(3).eq(0)
        return img.updateMask(mask).select('LST_Day_1km')
    else:
        qc = img.select('QC_Night')
        mask = qc.bitwiseAnd(3).eq(0)
        return img.updateMask(mask).select('LST_Night_1km')

def get_modis_for_date(target_date, roi_geom, region, out_folder, window=8):
    start = ee.Date((target_date - timedelta(days=window)).strftime('%Y-%m-%d'))
    end = ee.Date((target_date + timedelta(days=window)).strftime('%Y-%m-%d'))
    modis = ee.ImageCollection('MODIS/061/MOD11A1') \
        .filterDate(start, end) \
        .filterBounds(roi_geom)
    modis_list = modis.toList(modis.size())
    best_img = None
    min_diff = 9999
    
    # Skip if files already exist
    date_str = target_date.strftime('%Y-%m-%d')
    out_path_day = os.path.join(out_folder, f'MODIS_LST_Day_{date_str}.tif')
    out_path_night = os.path.join(out_folder, f'MODIS_LST_Night_{date_str}.tif')
    
    if os.path.exists(out_path_day) and os.path.exists(out_path_night):
        print(f"MODIS files for {date_str} already exist. Skipping.")
        # verify_image(out_path_day)
        # verify_image(out_path_night)
        return
    
    # Get size before the for loop to avoid repeated calls
    collection_size = modis.size().getInfo()
    if collection_size == 0:
        print(f"No MODIS data found for date {date_str}")
        return
        
    for i in range(collection_size):
        img = ee.Image(modis_list.get(i))
        img_date = datetime.strptime(img.date().format('YYYY-MM-dd').getInfo(), '%Y-%m-%d')
        diff = abs((img_date - target_date).days)
        if diff < min_diff:
            min_diff = diff
            best_img = img
    
    if best_img:
        lst_day = mask_modis_lst(best_img, day=True)
        lst_night = mask_modis_lst(best_img, day=False)
        
        # Use consistent CRS (EPSG:3857) for all exports
        export_ee_image(lst_day, ['LST_Day_1km'], region, out_path_day, 1000, crs='EPSG:3857')
        time.sleep(2)
        export_ee_image(lst_night, ['LST_Night_1km'], region, out_path_night, 1000, crs='EPSG:3857')
        time.sleep(2)

def main():
    base_dir = '/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_retrieval/new_download_data'
    export_base = '/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_retrieval/new_download_data'
    os.makedirs(export_base, exist_ok=True)
    
    for roi_folder in os.listdir(base_dir)[X:Y]:
        roi_path = os.path.join(base_dir, roi_folder)
        if os.path.isdir(roi_path):
            lst_dir = os.path.join(roi_path, 'lst')
            tif_files = list(Path(lst_dir).glob('*.tif'))
            if not tif_files:
                continue
                
            coords = get_roi_coords_from_tif(str(tif_files[0]))
            roi_geom = ee.Geometry.Polygon([coords], 'EPSG:3857')
            region = coords
            
            # Get original LST dates
            lst_dates = get_lst_dates(str(lst_dir))
            print(f"Found {len(lst_dates)} LST dates for {roi_folder}")
            
            # Generate missing dates between LST intervals
            missing_dates = generate_missing_dates(lst_dates)
            print(f"Identified {len(missing_dates)} missing dates to fill gaps")
            
            # Combine all dates to process
            all_dates = sorted(lst_dates + missing_dates)
            
            out_folder1 = os.path.join(export_base, roi_folder, 'era5')
            out_folder2 = os.path.join(export_base, roi_folder, 'modis')
            os.makedirs(out_folder1, exist_ok=True)
            os.makedirs(out_folder2, exist_ok=True)
            
            for date in all_dates:
                date_str = date.strftime('%Y-%m-%d')
                if date in lst_dates:
                    print(f'Processing LST date: {roi_folder}, Date: {date_str}')
                else:
                    print(f'Processing missing date: {roi_folder}, Date: {date_str} (gap-fill)')
                
                print(f'ERA5 for date: {date}')
                get_era5_for_date(date, roi_geom, region, out_folder1)
                
                print(f'MODIS for date: {date}')
                get_modis_for_date(date, roi_geom, region, out_folder2)
            
            # Check existing files for correct projection and fix if needed
            print(f"Verifying coordinate systems for {roi_folder}...")
            
            # Check ERA5 files
            era5_files = glob.glob(os.path.join(out_folder1, '*.tif'))
            for era5_file in era5_files:
                reproject_if_needed(era5_file, 'EPSG:3857')
                
            # Check MODIS files  
            modis_files = glob.glob(os.path.join(out_folder2, '*.tif'))
            for modis_file in modis_files:
                reproject_if_needed(modis_file, 'EPSG:3857')
            
            print(f"Finished processing {roi_folder}")
            print("--------------------------------")

if __name__ == '__main__':
    main()
