#%%
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.enums import Resampling
import os
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import geopandas as gpd

import parameters as param


class download_scene():
    '''
    This is an empty class to mention that the images could be dowloaded dynamically. 
    This sectin could be develop in a later time.

    '''
    def __init__(self):
        pass



class Raster_processing():
    # Dictionnary which contains the necessary bands to perform the requested spectral index. The dictionnary can contain
    # several indices such as NDVI, NDWI, etc.
    bands_to_select_dict = {'S2_NBR': ('B02_10m','B03_10m','B04_10m','B08_10m','B12_20m','SCL_20m'),
                            'S2_NDVI': ('B02_10m','B03_10m','B04_10m','B08_10m'),
                            }


    def __init__(self, raster_path, sensor_name: str, spectral_index_name: str):
        self.raster_to_read = raster_path
        self.sensor = sensor_name
        self.spectral_index = spectral_index_name


    def open_raster(self, raster_path):
        '''
        Open rasterio object to prepare for the reading. The goal is to open all raster first, then read them to match the specific needs.
        Ex.: reading multiple scenes and their respective bands and wanting to match the same bands at the moment of reading.
        '''
        raster_open_object = []     

        for path in raster_path:

            raster_open_object.append(rasterio.open(path, mode="r"))

        return raster_open_object

    def read_raster(self):

        pass

    def get_band_path(self, bands_to_use):

        band_path_process_selection = {'S2': self.parse_manifest_xml(bands=bands_to_use),}

        return band_path_process_selection[self.sensor]

    def parse_manifest_xml(self, bands, manifest_path=None):

        '''
        Parse the Manifest.safe to get the path of Sentinel-2 bands within the folder structure of the scene.

        manifest_path: 
        bands: Band names used to query the manifest.safe (xml) file.

        '''

        if not manifest_path:
            manifest_path = os.path.join(self.raster_to_read,'manifest.safe')


        self.selected_bands_path = []

        tree = ET.parse(manifest_path)
        root = tree.getroot()

        for element in root.iter('dataObject'):

            # print(element.attrib['ID'])

            if any(e in element.attrib['ID'] for e in bands):                
                for child in element:
                    self.selected_bands_path.append(child.find('fileLocation').attrib['href'])

        return self.selected_bands_path
    
    def select_bands(self):

        # Build the dictionnary key to select the necessary bands
        selection_key_spectral_index = self.sensor + '_' + self.spectral_index

        self.bands_to_select = self.bands_to_select_dict[selection_key_spectral_index]

        return self.bands_to_select

#read raster
#resampling
#export raster

class index_calculation:

    def __init__(self):
        pass

#NBR
#NDVI
#index threshold as variable not function

class Raster_IO:
    
    

    def __init__(self, folder_path):

        self.folder_path = folder_path
        self.list_file_scene = os.listdir(self.folder_path)
        self.abs_path_list_file_scene = self.create_absolute_path()
        self.raster_path_to_process = []

# insert error raise if no folder is provided.

    def extract_raster(self):

        self.check_if_file_zip()

        #return self.raster_path_to_process

    def get_raster_paths(self):
        return self.raster_path_to_process

    def check_if_file_zip(self):
        '''
        Validates if one or several .zip are located within the repository. If a .zip file is detected, it is automatically unzipped.
        In the case the .zip have already been unzipped, then it won't be unzipped.
        '''
        if len(self.abs_path_list_file_scene) == 1 and os.path.splitext(self.abs_path_list_file_scene[0])[1] =='.zip':

            
            self.unzip(self.abs_path_list_file_scene[0])
            os.path.splitext(self.raster_path_to_process.append(self.abs_path_list_file_scene[0]))
        
        elif len(self.abs_path_list_file_scene) == 1 and os.path.splitext(self.abs_path_list_file_scene[0])[1] =='.SAFE':
            
             self.raster_path_to_process.append(self.abs_path_list_file_scene[0])


        elif len(self.abs_path_list_file_scene) > 1:
            
            list_file_scene_roots = [os.path.splitext(file)[0] for file in self.abs_path_list_file_scene if os.path.splitext(file)[1] == '.zip']
            
            for root in list_file_scene_roots:

                if not os.path.exists(root):

                    self.unzip(root+'.zip')
                    self.raster_path_to_process.append(root)
     
                else:
                    self.raster_path_to_process.append(root)
        else:
            raise FileNotFoundError("File not recognized or doesn't exist. Please verify if a .zip or .SAFE file(s) are in the folder.")
            


    def create_absolute_path(self):
        paths = []
        for file in self.list_file_scene:
            paths.append(os.path.join(self.folder_path, file))

        return paths
    

    def unzip(self, file_to_unzip):

        with ZipFile(file_to_unzip, 'r') as unzip_file:
            unzip_file.extractall(path=self.folder_path)
        
def get_processing_bound(path):
    '''
    Calculate the total bounds of a geometry.

    path: Path to geometry file. For the moment, the crs should be the same as the input rasters. In later improvements, the crs will automatically be managed

    Return: Tuple containing (minX, minY, maxX, maxY)
    '''

    geom = gpd.read_file(path)

    bound_geom = geom.total_bounds
    
    return tuple(bound_geom)


def create_mask(mask_source: np.array, target_value):
    '''
    Create a boolean mask to mask out the water bodies.

    mask_source: array which is used to create the mask.

    target_value: value(s) to be used to create the mask.

    return boolean mask
    '''

    # mask = np.where(mask_source != mask_value, 1, 0).astype(np.uint8)
    
    
    mask = mask_source != target_value

    return mask
    


def apply_mask(array_to_mask, mask_array, nodata):
    
    return array_to_mask * mask_array + ~mask_array * nodata

    

def calculate_S2_NBR(band_08, band_12):
    band_08_scale = band_08/10000 
    band_12_scale = band_12/10000

    NBR = (band_08_scale - band_12_scale)/(band_08_scale + band_12_scale)
    
    del band_08_scale, band_12_scale

    return NBR

def reclassify_NBR():

    pass


#%%

def main():
        
    scene_repositary = param.scene_rep[0]

    geom_file = param.geom_file[0]
    sensor = param.sensor
    spectral_index = param.spectral_index
    no_data = param.no_data


    # Dict to indicate the order of bands that are used for specific calculated indices.
    # When an indice is selected, the sub dictionnary is used as a guide to stack the bands to be used.
    band_order_for_indices = {'S2_NBR':('B02','B03', 'B04','B08','B12','SCL'),
                              'S2_NDVI': ('B02','B03', 'B04','B08')
                                }  

    process_name = sensor + '_' + spectral_index


    for key, repositary in scene_repositary.items():
        
        scene_name_file = Raster_IO(repositary)
        scene_name_file.extract_raster()

        print(repositary)
        
        open_raster_list = []

        # Process each raster scene to get band path, open the rasters, merge tiles, clip, resample, calculate indices.
        for path in scene_name_file.get_raster_paths():
            print(path, '\n')
            scene_processing = Raster_processing(path, sensor_name=sensor, spectral_index_name=spectral_index)

            bands_to_process = scene_processing.select_bands()
            band_path = scene_processing.get_band_path(bands_to_use=bands_to_process)
            band_path = [os.path.join(path, p) for p in band_path]
           
            open_raster = scene_processing.open_raster(raster_path=band_path)
            open_raster_list.append(open_raster)
            del scene_processing, bands_to_process, band_path, open_raster
            # Here I assume that all rasters have the same coordinate system. An additional step to verify and test if they all have the same.
            # For example, in case of difference, we could assign to all EPSG:4326. This would make the process more robust. 

            print(open_raster_list)

        #del scene_name_file
        # Get bounds of processing area
        processing_bounds = get_processing_bound(path=geom_file)

        # Dict to store extracted bands for further process
        dict_keys = band_order_for_indices[process_name]
        extracted_bands = dict.fromkeys(dict_keys)
        
        # Pairing band pairs for each scenes, in the goal of performing further processing on merged scenes band per band.
        for band in range(0, len(open_raster_list[0])):
            
            raster_to_merge = []

            for index in range(0, len(open_raster_list)):

                raster_to_merge.append(open_raster_list[index][band])
                
            # bounds is used to limit the merge to the total extent of the provided geometry.
            # This helps improving the rapidity of processing since non-necessary data is not processed.

            if process_name == 'S2_NBR' and dict_keys[band] == 'B12':

                band_merged, out_transform = merge(raster_to_merge, bounds= processing_bounds, res=10, resampling=Resampling.bilinear)

                extracted_bands[dict_keys[band]] = band_merged


            elif process_name == 'S2_NBR' and dict_keys[band] == 'SCL':
                
                band_merged, out_transform = merge(raster_to_merge, bounds= processing_bounds, res=10, resampling=Resampling.nearest)

                extracted_bands[dict_keys[band]] = band_merged


            elif process_name == 'S2_NBR':

                band_merged, out_transform = merge(raster_to_merge, bounds= processing_bounds, resampling=Resampling.bilinear)

                extracted_bands[dict_keys[band]] = band_merged


            del raster_to_merge, band_merged


        # Create boolean mask to exclude areas such as water bodies. It would also be possible to include other classes.
        mask = create_mask(mask_source=extracted_bands['SCL'], target_value=6)
        
        # Calculating the NBR from S2 B09 and B12 bands.
        index_NBR = calculate_S2_NBR(extracted_bands['B08'], extracted_bands['B12'])
        
        # Converting the masked pixels to nodata
        index_NBR_mask = apply_mask(index_NBR, mask, nodata=no_data)
        
        del extracted_bands

        # Preparing to export the NBR raster if it is wished to see the results.
        output_meta = open_raster_list[0][0].meta.copy()
        output_meta.update({"driver": "GTIFF",
                            "height": index_NBR_mask.shape[1],
                            "width": index_NBR_mask.shape[2],
                            "count": 1,
                            "transform": out_transform,
                            "dtype": 'float64',
                            "nodata": no_data             
                            })
        
        
        
        with rasterio.open(os.path.join(repositary, f"{key}_NBR_index.tif"), mode='w', **output_meta,) as dst:

            dst.write(index_NBR_mask[0], 1)

        del open_raster_list, mask, index_NBR, index_NBR_mask, #scene_name_file
 
if __name__ == '__main__':

    main()



# # Enables the option to clip the reading window to boundaries from a geometry file.
#        
# def open_raster(self, raster_path, geom_file_path: str = None, clip_to_window: bool = False):
#  if clip_to_window:
#             read_window_bounds = self.get_processing_bound(geom_file_path)

#             with rasterio.open(pat[h, mode="r", window=rasterio.windows.from_bounds(read_window_bounds)) as src:

# %%
