# import os
# import rasterio
# import shutil # Added for file operations if needed by reproject_if_needed
# from datetime import datetime, timedelta # Added for potential use if adapting functions
# # Assuming new_main.py is in the same directory or Python path
# from new_main import reproject_if_needed, verify_image 

# path = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_retrieval/download_data"
# location = "LeLoi_KienXuong_ThaiBinh" # Example location

# era5_path = os.path.join(path, location, "era5/ERA5_temperature_skin_2022-12-20.tif")
# modis_path = os.path.join(path, location, "modis/MODIS_LST_Day_2022-12-20.tif")
# ndvi_path = os.path.join(path, location, "LeLoi_ndvi8days/ndvi8days_2022-12-23.tif")

# # target_crs = 'EPSG:3857'

# # print("--- Initial State ---")
# # print("ERA5:")
# # verify_image(era5_path)
# # print("\nMODIS:")
# # verify_image(modis_path)
# # print("\nNDVI:")
# # verify_image(ndvi_path)

# # print(f"\n--- Attempting to reproject NDVI to {target_crs} if needed ---")
# # # reprojected = reproject_if_needed(ndvi_path, target_crs)

# # if reprojected:
# #     print(f"NDVI file was reprojected to {target_crs}\n")
# # else:
# #     print(f"NDVI file already has CRS {target_crs} or an error occurred during reprojection check.\n")

# # print("--- Final State ---")
# # print("ERA5:")
# # verify_image(era5_path)
# # print("\nMODIS:")
# # verify_image(modis_path)
# # print("\nNDVI (potentially reprojected):")
# # verify_image(ndvi_path)


# # Original code commented out for clarity, replaced by verify_image calls
# # Open the TIF file
# with rasterio.open(era5_path) as src:
#     # Read the image data
#     image_data = src.read() # Reading data not necessary for verification
    
#     # Get metadata
#     print(f"Image dimensions: {src.shape}")
#     print(f"Number of bands: {src.count}")
#     print(f"Coordinate system: {src.crs}")
#     print(f"Transform: {src.transform}")

# with rasterio.open(modis_path) as src:
#     # Read the image data
#     image_data = src.read()
    
#     # Get metadata
#     print(f"Image dimensions: {src.shape}")
#     print(f"Number of bands: {src.count}")
#     print(f"Coordinate system: {src.crs}")
#     print(f"Transform: {src.transform}")
    
# with rasterio.open(ndvi_path) as src: # This would read the potentially reprojected file
#     # Read the image data
#     image_data = src.read()
    
#     # Get metadata
#     print(f"Image dimensions: {src.shape}")
#     print(f"Number of bands: {src.count}")
#     print(f"Coordinate system: {src.crs}")
#     print(f"Transform: {src.transform}")



import os


path = "/mnt/hdd12tb/code/nhatvm/BRIOS/BRIOS/data_crawled"

for i in range(len(list(os.listdir(path)))):
    if list(os.listdir(path))[i] == "HoaKhanh_BuonMaThuot_DakLak":
        print(i)

