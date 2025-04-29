import ee, math, os
import requests
import shutil
import time
import zipfile
import rasterio
import numpy as np
from rasterio.merge import merge

ee.Initialize(project='ee-hadat')

# --- Helper Functions ---

def has_data_in_roi(image, roi):
    bands = ['VV', 'VH', 'RVI']
    for band in bands:
        if band not in image.bandNames().getInfo():
            return False
        stats = image.select(band).reduceRegion(ee.Reducer.minMax(), geometry=roi, scale=10)
        min_val = stats.get(band + '_min')
        if min_val is None:
            return False
    return True

def edge_correction(image):
    """
    Applies an edge correction to remove pixels with no data along the image boundary.
    """
    edge_mask = image.mask().reduce(ee.Reducer.allNonZero())
    return image.updateMask(edge_mask).copyProperties(image, ['system:time_start'])

def get_sentinel1_collection(start_date, end_date, roi):
    """
    Loads and clips the Sentinel-1 GRD collection filtered by date, region, and polarization.
    """
    s_date = start_date.advance(-8, 'day')
    e_date = end_date.advance(8, 'day')
    
    collection = (ee.ImageCollection('COPERNICUS/S1_GRD')
                  .filterBounds(roi)
                  .filterDate(s_date, e_date)
                  .filter(ee.Filter.eq('instrumentMode', 'IW'))
                  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')))
    
    def process(image):
        image = image.clip(roi)
        image = edge_correction(image)
        return image

    return collection.map(process)

def calculate_rvi(image):
    """
    Calculates the Radar Vegetation Index (RVI) from the VV and VH bands.
    Formula: RVI = (4 * VH) / (VV + VH)
    """
    vv = image.select('VV')
    vh = image.select('VH')
    rvi = vh.multiply(4).divide(vv.add(vh)).rename('RVI')
    return image.addBands(rvi).set('system:time_start', image.get('system:time_start'))

def calculate_8day_composites_sar(image_collection, start_date, end_date):
    days_step = 8
    start = ee.Date(start_date)
    end = ee.Date(end_date)
    millis_step = days_step * 24 * 60 * 60 * 1000
    list_of_dates = ee.List.sequence(start.millis(), end.millis(), millis_step)

    def composite_for_millis(millis):
        composite_center = ee.Date(millis)
        composite_start = composite_center.advance(-8, 'day')
        composite_end = composite_center.advance(8, 'day')
        period_collection = image_collection.filterDate(composite_start, composite_end)
        composite = ee.Image(ee.Algorithms.If(
            period_collection.size().gt(0),
            calculate_rvi(period_collection.median()).set('system:time_start', composite_center.millis()),
            ee.Image(0).updateMask(ee.Image(0)).set('system:time_start', composite_center.millis())
        ))
        return composite

    composites = ee.ImageCollection(list_of_dates.map(composite_for_millis))
    return composites

def sort_by_time(composites):
    return composites.sort('system:time_start')

def smooth_time_series(composites):
    image_list = composites.toList(composites.size())
    collection_size = composites.size().getInfo()

    def has_rvi(img):
        return ee.Number(img.bandNames().size()).gt(3)

    smoothed_images = []
    for i in range(collection_size):
        image = ee.Image(image_list.get(i))
        image_date = ee.Date(image.get('system:time_start'))
        previous = ee.Image(image_list.get(i - 1)) if i > 0 else image
        next_img = ee.Image(image_list.get(i + 1)) if i < (collection_size - 1) else image

        current_has = has_rvi(image)
        previous_has = has_rvi(previous)
        next_has = has_rvi(next_img)

        smoothed = ee.Image(ee.Algorithms.If(
            current_has,
            ee.Image(ee.Algorithms.If(
                previous_has.And(next_has),
                ee.ImageCollection([previous, image, next_img]).mean().set('system:time_start', image_date.millis()),
                ee.Image(ee.Algorithms.If(
                    next_has,
                    ee.ImageCollection([image, next_img]).mean().set('system:time_start', image_date.millis()),
                    ee.Image(ee.Algorithms.If(
                        previous_has,
                        ee.ImageCollection([image, previous]).mean().set('system:time_start', image_date.millis()),
                        image
                    ))
                ))
            )),
            ee.Image(ee.Algorithms.If(
                previous_has.And(next_has),
                ee.ImageCollection([previous, next_img]).mean().set('system:time_start', image_date.millis()),
                ee.Image(ee.Algorithms.If(
                    previous_has,
                    previous.set('system:time_start', image_date.millis()),
                    ee.Image(ee.Algorithms.If(
                        next_has,
                        next_img.set('system:time_start', image_date.millis()),
                        image
                    ))
                ))
            ))
        ))
        smoothed_images.append(smoothed)

    return ee.ImageCollection(smoothed_images)

def download_band(image, band_name, roi, date_str, out_folder, max_retries=4):
    local_path = os.path.join(out_folder, f"{band_name}_{date_str}.tif")
    temp_dir = os.path.join(out_folder, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    for attempt in range(max_retries):
        try:
            url = image.select([band_name]).getDownloadURL({
                'scale': 10,
                'region': roi,
                'fileFormat': 'GeoTIFF',
                'maxPixels': 1e13,
                'expires': 3600
            })
            print(f"Attempt {attempt+1} for {band_name} {date_str}: {url}")

            temp_zip_path = os.path.join(temp_dir, "download.zip")
            response = requests.get(url, stream=True, timeout=300, headers={'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64)'})
            print(f"Status code: {response.status_code}")
            print(f"Response headers: {response.headers}")

            response.raise_for_status()

            with open(temp_zip_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        file.write(chunk)

            with zipfile.ZipFile(temp_zip_path, 'r') as zip_ref:
                print(f"Zip file contents: {zip_ref.namelist()}")
                zip_ref.extractall(temp_dir)

            print(f"Files in temp dir: {os.listdir(temp_dir)}")
            tif_files = [f for f in os.listdir(temp_dir) if f.endswith('.tif')]

            if len(tif_files) == 1:
                tif_file = tif_files[0]
                src_path = os.path.join(temp_dir, tif_file)
                shutil.copy(src_path, local_path)  # Use shutil.copy instead of os.replace to handle cross-device links
                if is_valid_tif(local_path):
                    print(f"Successfully downloaded {band_name} for {date_str}")
                    # Clean up temp directory after successful download
                    shutil.rmtree(temp_dir)
                    return local_path
                else:
                    print(f"⚠️ Invalid file for {band_name} {date_str}, retrying...")
                    time.sleep(2 * (attempt + 1))
            else:
                print(f"⚠️ Unexpected number of .tif files in zip: {len(tif_files)}")
                time.sleep(2 * (attempt + 1))
        except requests.exceptions.RequestException as e:
            print(f"⚠️ Request error for {band_name} {date_str}: {e}")
            time.sleep(2 * (attempt + 1))
        except Exception as e:
            print(f"Error downloading {band_name} for {date_str}: {e}")
            time.sleep(2 * (attempt + 1))
    
    print(f"Failed to download {band_name} for {date_str} after {max_retries} attempts.")
    # Clean up temp directory on failure
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    return None

def merge_bands(rvi_path, vv_path, vh_path, output_path):
    """Merge RVI, VV, VH into a single 3-band image (reshaping if necessary)."""
    band_paths = [rvi_path, vv_path, vh_path]
    bands = []

    for path in band_paths:
        with rasterio.open(path) as src:
            band = src.read(1)
            bands.append((band, src.profile))

    min_shape = min(band[0].shape for band in bands)
    bands_resized = [np.resize(band[0], min_shape) for band in bands]

    profile = bands[0][1]
    profile.update(count=3, dtype='float32')

    with rasterio.open(output_path, "w", **profile) as dst:
        for i, band in enumerate(bands_resized, start=1):
            dst.write(band, i)

    print(f"✅ Merged 3-band image saved: {output_path}")

def is_valid_tif(file_path):
    if not os.path.exists(file_path) or os.path.getsize(file_path) < 1024:
        print(f"Invalid: {file_path} - File missing or too small")
        return False
    try:
        with rasterio.open(file_path) as src:
            if src.count == 0 or src.width == 0 or src.height == 0:
                print(f"Invalid: {file_path} - No bands or zero dimensions")
                return False
            return True
    except rasterio.errors.RasterioIOError as e:
        print(f"Invalid: {file_path} - Rasterio error: {e}")
        return False

def export_sentinel1_rvi(sentinel_collection, big_folder, roi, image_name, roi_name, folder_name):
    out_folder = os.path.join(big_folder, roi_name, folder_name)
    os.makedirs(out_folder, exist_ok=True)

    image_list = sentinel_collection.toList(sentinel_collection.size())
    count = sentinel_collection.size().getInfo()

    for i in range(count):
        image = ee.Image(image_list.get(i))
        date_str = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd').getInfo()
        
        if not has_data_in_roi(image, roi):
            print(f"No data for {date_str}, skipping.")
            continue

        rvi_path = download_band(image, 'RVI', roi, date_str, out_folder)
        time.sleep(2)
        vv_path = download_band(image, 'VV', roi, date_str, out_folder)
        time.sleep(2)
        vh_path = download_band(image, 'VH', roi, date_str, out_folder)
        time.sleep(3)

        if rvi_path and vv_path and vh_path:
            output_path = os.path.join(out_folder, f"{image_name}_{date_str}.tif")
            merge_bands(rvi_path, vv_path, vh_path, output_path)
            # Remove individual band files after successful merge
            for band_path in [rvi_path, vv_path, vh_path]:
                if os.path.exists(band_path):
                    os.remove(band_path)
                    print(f"Removed individual band file: {band_path}")
        else:
            print(f"⚠️ Skipping merge for {date_str} due to missing or invalid bands.")
            # Clean up any downloaded files if merge fails
            for band_path in [rvi_path, vv_path, vh_path]:
                if band_path and os.path.exists(band_path):
                    os.remove(band_path)
                    print(f"Removed partial band file: {band_path}")

def display_rvi(rvi_collection, layer_name):
    image_list = rvi_collection.toList(rvi_collection.size())
    count = rvi_collection.size().getInfo()
    for i in range(count):
        rvi_image = ee.Image(image_list.get(i))
        if rvi_image.bandNames().contains('RVI').getInfo():
            date = ee.Date(rvi_image.get('system:time_start')).format('YYYY-MM-dd').getInfo()
            print(f"{layer_name} - Image with RVI from {date}")

def export_sentinel1_rvi_drive(sentinel_collection, roi, image_name, folder_name):
    image_list = sentinel_collection.toList(sentinel_collection.size())
    count = sentinel_collection.size().getInfo()
    for i in range(count):
        image = ee.Image(image_list.get(i))
        date_str = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd').getInfo()
        task = ee.batch.Export.image.toDrive(
            image=image.select(['RVI', 'VV', 'VH']),
            description=image_name + date_str,
            folder=folder_name,
            scale=10,
            region=roi,
            fileFormat='GeoTIFF',
            maxPixels=1e13
        )
        task.start()
        print('Export initiated for image on date:', date_str)


def main_rvi(start_date, end_date, ROI, big_folder, roi_name):
    sentinel1_collection = get_sentinel1_collection(start_date, end_date, ROI)
    print('Sentinel-1 collection size:', sentinel1_collection.size().getInfo())

    sentinel_8day_composites = calculate_8day_composites_sar(sentinel1_collection, start_date, end_date)
    print('8-day composites size:', sentinel_8day_composites.size().getInfo())

    sorted_composites = sort_by_time(sentinel_8day_composites)
    print('First composite band names:', sorted_composites.first().bandNames().getInfo())

    smoothed_composites = smooth_time_series(sorted_composites)

    export_sentinel1_rvi(smoothed_composites, big_folder, ROI, 'rvi_8days', roi_name, f'{roi_name.split("_")[0]}_rvi_8days')
    # export_sentinel1_rvi_drive(smoothed_composites, ROI, 'rvi_8days_', 'LST')

# # --- Main Execution ---

# start_date = ee.Date('2023-01-01')
# end_date = ee.Date('2023-02-01')
# ROI = ee.Geometry.Polygon([
#     [[11855180.0, 2328782.140517925], [11855180.0, 2335016.8225009246], [11851160.0, 2335016.8225009246], [11851160.0, 2328782.140517925], [11855180.0, 2328782.140517925]]
# ], 'EPSG:3857')

# big_folder = "/mnt/data1tb/LSTRetrieval/Code/download_data"

# sentinel1_collection = get_sentinel1_collection(start_date, end_date, ROI)
# print('Sentinel-1 collection size:', sentinel1_collection.size().getInfo())

# sentinel_8day_composites = calculate_8day_composites_sar(sentinel1_collection, start_date, end_date)
# print('8-day composites size:', sentinel_8day_composites.size().getInfo())

# sorted_composites = sort_by_time(sentinel_8day_composites)
# print('First composite band names:', sorted_composites.first().bandNames().getInfo())

# smoothed_composites = smooth_time_series(sorted_composites)

# export_sentinel1_rvi(smoothed_composites, big_folder, ROI, 'rvi_8days', "BinhThanh_DucHue_LongAn", 'BinhThanh_rvi_8days')
# # export_sentinel1_rvi_drive(smoothed_composites, ROI, 'rvi_8days_', 'LST')