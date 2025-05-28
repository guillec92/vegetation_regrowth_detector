# Vegetation Regrowth Detector After Wildfire

This processor analyses Sentinel-2 (L2A) scenes to detect burnt areas and assess the intensity of the vegetation regrowth.

The processor takes as input 3 scenes (prefire, postfire (immediately after) and a recent scene). The burnt areas and vegetation regrowth are measured by using the Normalized Burnt Ratio (NBR) indice.

## Input parameters

The folder path of the scenes must be in full path.

ex.: "/home/user/../20240910"

Under the "20240910" folder, the scenes are stored.

## Output results

* before_fire_S2_NBR.tif
  * Value of NBR before the fire

* after_fire_S2_NBR.tif
  * Value of NBR after the fire
    
* recent_S2_NBR.tif
  * Value of NBR in the recent days after the fire


* burnt_area.tif
  * Calculation of dNBR (prefire NBR - postfire NBR)
    
    1: Low-severity burn
    
    2: Moderate-low severity burn
    
    3: Moderate-high severity burn
    
    4: High-severity burn
    
    255: nodata
    

* regrowth_masked.tif
  * Calculation of dNBR (postfire NBR - recent NBR)

    1: High post-fire regrowth
    
    2: Low post-fire regrowth
    
    3: No regrowth
    
    255: nodata
    
