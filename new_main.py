import os
import glob
import time
from datetime import datetime, timedelta
from pathlib import Path
import rasterio
from rasterio.warp import transform_bounds
import ee
from rvi_retrieval import get_sentinel1_collection, calculate_8day_composites_sar, sort_by_time, smooth_time_series

# Initialize Earth Engine
try:
    ee.Initialize()
except Exception:
    ee.Authenticate()
    ee.Initialize()

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

def find_nearest_date(target_date, date_list, window=8):
    min_diff = timedelta(days=window+1)
    nearest = None
    for d in date_list:
        diff = abs((d - target_date).days)
        if diff <= window and diff < min_diff.days:
            min_diff = timedelta(days=diff)
            nearest = d
    return nearest

def export_ee_image(image, bands, region, out_path, scale):
    # Create proper geometry object with explicit projection
    region_geometry = ee.Geometry.Polygon([region], 'EPSG:3857')
    
    # Clip the image first and make sure it has data
    image = image.clip(region_geometry)
    
    # Then create download URL with proper parameters
    url = image.select(bands).getDownloadURL({
        'scale': scale,
        'region': region_geometry,
        'fileFormat': 'GeoTIFF',
        'maxPixels': 1e13,
        'format': 'GEO_TIFF',
        'expires': 3600,
    })
    
    import requests, zipfile, shutil
    temp_dir = os.path.dirname(out_path) + '/temp_dl'
    os.makedirs(temp_dir, exist_ok=True)
    temp_zip_path = os.path.join(temp_dir, 'download.zip')
    try:
        r = requests.get(url, stream=True, timeout=300, headers={'User-Agent': 'Mozilla/5.0'})
        # Check if the response is valid before saving
        if r.status_code != 200:
            raise Exception(f"Bad response: {r.status_code}, {r.text}")
            
        with open(temp_zip_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    
        # Verify file is a zip before attempting to extract
        if zipfile.is_zipfile(temp_zip_path):
            with zipfile.ZipFile(temp_zip_path, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            tif_files = [f for f in os.listdir(temp_dir) if f.endswith('.tif')]
            if tif_files:
                shutil.move(os.path.join(temp_dir, tif_files[0]), out_path)
        else:
            print(f"Downloaded file is not a valid ZIP: {temp_zip_path}")
            with open(temp_zip_path, 'rb') as f:
                print(f"First 100 bytes: {f.read(100)}")
    except Exception as e:
        print(f"Error downloading/exporting {out_path}: {e}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def get_rvi_for_date(target_date, roi_geom, region, out_folder, window=8):
    start = ee.Date((target_date - timedelta(days=window)).strftime('%Y-%m-%d'))
    end = ee.Date((target_date + timedelta(days=window)).strftime('%Y-%m-%d'))
    collection = get_sentinel1_collection(start, end, roi_geom)
    composites = calculate_8day_composites_sar(collection, start, end)
    sorted_composites = sort_by_time(composites)
    smoothed = smooth_time_series(sorted_composites)
    image_list = smoothed.toList(smoothed.size())
    best_img = None
    min_diff = 9999
    for i in range(smoothed.size().getInfo()):
        img = ee.Image(image_list.get(i))
        img_date = datetime.utcfromtimestamp(img.get('system:time_start').getInfo()/1000)
        diff = abs((img_date - target_date).days)
        if diff < min_diff:
            min_diff = diff
            best_img = img
    if best_img:
        out_path = os.path.join(out_folder, 'RVI.tif')
        export_ee_image(best_img, ['RVI', 'VV', 'VH'], region, out_path, 10)
        time.sleep(2)

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
    out_path = os.path.join(out_folder, 'ERA5_temperature_skin.tif')
    export_ee_image(best_img, ['temperature_2m', 'skin_temperature'], region, out_path, 11132)
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
    for i in range(modis.size().getInfo()):
        img = ee.Image(modis_list.get(i))
        img_date = datetime.strptime(img.date().format('YYYY-MM-dd').getInfo(), '%Y-%m-%d')
        diff = abs((img_date - target_date).days)
        if diff < min_diff:
            min_diff = diff
            best_img = img
    if best_img:
        lst_day = mask_modis_lst(best_img, day=True)
        lst_night = mask_modis_lst(best_img, day=False)
        out_path_day = os.path.join(out_folder, 'MODIS_LST_Day.tif')
        out_path_night = os.path.join(out_folder, 'MODIS_LST_Night.tif')
        export_ee_image(lst_day, ['LST_Day_1km'], region, out_path_day, 1000)
        time.sleep(2)
        export_ee_image(lst_night, ['LST_Night_1km'], region, out_path_night, 1000)
        time.sleep(2)

def main():
    base_dir = 'download_data'
    export_base = 'exported_data'
    os.makedirs(export_base, exist_ok=True)
    for roi_folder in os.listdir(base_dir):
        roi_path = os.path.join(base_dir, roi_folder)
        if os.path.isdir(roi_path):
            lst_dir = os.path.join(roi_path, 'lst')
            tif_files = list(Path(lst_dir).glob('*.tif'))
            if not tif_files:
                continue
            coords = get_roi_coords_from_tif(str(tif_files[0]))
            roi_geom = ee.Geometry.Polygon([coords], 'EPSG:3857')
            region = coords
            lst_dates = get_lst_dates(str(lst_dir))
            for date in lst_dates:
                date_str = date.strftime('%Y-%m-%d')
                out_folder = os.path.join(export_base, roi_folder, date_str)
                os.makedirs(out_folder, exist_ok=True)
                print(f'Exporting for ROI: {roi_folder}, Date: {date_str}')
                print(f'RVI for date: {date}')
                get_rvi_for_date(date, roi_geom, region, out_folder)
                print(f'ERA5 for date: {date}')
                get_era5_for_date(date, roi_geom, region, out_folder)
                print(f'MODIS for date: {date}')
                get_modis_for_date(date, roi_geom, region, out_folder)

if __name__ == '__main__':
    main()
