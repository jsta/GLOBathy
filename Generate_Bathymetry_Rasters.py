# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Name:         Generate_Bathymetry_Rasters
# Purpose:      This script is designed to calculate bathymetric maps for a waterbodies in Tagged 
#               Image File Format (TIFF) and in WGS84 projection system based on a simple equation
#               of linear distance from shoreline. Two inputs are required:
#                   1) a csv file that includes the maximum depths of the waterbodies
#                   e.g., “GLOBathy_basic_parameters(ALL_LAKES).csv” can be used as a template
#                   2) polygon shapefiles of the corresponding waterbodies
#                   e.g., “HydroLAKES_polys_v10.shp” of the HydroLAKES dataset can be used as a template
#               This script can be used to generate bathymetry rasters for a waterbody that is not
#               available in the GLOBathy dataset or if updated observation/estimate of maximum depth
#               for a waterbody becomes available.
# Created:      2020
#
# Notes:        The open-source GDAL python library was used in this script to perform the distance
#               calculations using the "ComputeProximity" function against a gridded lake mask:
#                   GDAL/OGR contributors (2020). GDAL/OGR Geospatial Data Abstraction software 
#                   Library. Open Source Geospatial Foundation. URL https://gdal.org
#
# Disclaimer:   THIS CODE IS PROVIDED “AS IS”. IN NO EVENT SHALL PAGERDUTY OR CONTRIBUTORS BE LIABLE
#               FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, ORCONSEQUENTIAL DAMAGES
#               (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
#               USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) SUSTAINED BY YOU OR A THIRD PARTY,
#               HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#               OR TORT ARISING IN ANY WAY OUT OF THE USE OF THIS CODE, EVEN IF ADVISED OF THE
#               POSSIBILITY OF SUCH DAMAGE.
#
# Citation:     Khazaei, B., Read, L. K., Casali, M., Sampson, K. M. & Yates, D. N.
#               GLOBathy, the Global Lakes Bathymetry Dataset. figshare (2021).
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

# Import Python core modules
import os
import sys
import time

# Import Additional Modules
import ogr
import osr
import gdal
from gdalconst import *
import numpy
import pandas as pd

# Assists in using BandWriteArray, BandReadAsArray, and CopyDatasetInfo
from gdalnumeric import *

sys.dont_write_bytecode = True
print('Script initiated at {0}'.format(time.ctime()))

# Module configurations
gdal.UseExceptions()  # this allows GDAL to throw Python Exceptions
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.AllRegister()

# Input lake vector layer, driver and field name information
inLakes = "/path/to/file.shp"
inDriverName = 'ESRI Shapefile'         # 'GPKG' / 'ESRI Shapefile'
fieldname = 'Hylak_id'                  # The fieldname to identify a lake polygon

# Input lake depth information
lake_depths = "/path/to/file.csv"       # Maximum depth CSV
lake_depths_maxdepth_field = 'Dmax_use_m'  # CSV column containing maximum depth
lake_depths_id_field = 'Hylak_id'       # CSV column containing lake ID

# Output directory, format, suffix, and nodata information
out_tiles_dir = "path/to/output/dir"
envelopeSHP = os.path.join(out_tiles_dir, 'Raster_Tiles.shp')  # Output raster tile boundary shapefile (optional)
outFormat = 'GTiff'                     # Control the output format driver for raster tiles
suffix = 'tif'                          # Filename suffix to use for output filenames
NoDataVal = -9999                       # Nodata value to give all non-lake depth pixels

# Define the output grid CRS, bounds, and cell size
outGridProj = 4326
demProj = osr.SpatialReference()        # Instantiate the spatial reference object
demProj.ImportFromEPSG(outGridProj)     # Define using EPSG code specified in header
cellsize = float(1/3600)                # 1 arc-second expressed as a fraction of 1 degree.
x_min = -180.0                          # Origin: Lower Left X
y_min = -90.0                           # Origin: Lower Left Y
x_res = int(360.0/cellsize)
y_res = int(180.0/cellsize)

# Options for creating output files
writeenvelopeSHP = False                # Option to produce a shapefile of the subsetted DEM boundaries
copyIntermediateToDisk = False          # Save the rasterized polygon to disk?

# ---------- Functions ---------- #
def checkfield(layer, fieldname, string1):
    """
    Check for existence of provided fieldnames
    """
    layerDefinition = layer.GetLayerDefn()
    fieldslist = []
    for i in range(layerDefinition.GetFieldCount()):
        fieldslist.append(layerDefinition.GetFieldDefn(i).GetName())
    if fieldname in fieldslist:
        i = fieldslist.index(fieldname)
        field_defn = layerDefinition.GetFieldDefn(i)
    else:
        print('    Field {0} not found in input {1}. Terminating...'.format(fieldname, string1))
        raise SystemExit
    return field_defn, fieldslist

def my_envelope(envelope):
    """
    Creates a polygon object from an OGR envelope.
    """
    (minX, maxX, minY, maxY) = envelope

    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(minX, minY)
    ring.AddPoint(maxX, minY)
    ring.AddPoint(maxX, maxY)
    ring.AddPoint(minX, maxY)
    ring.AddPoint(minX, minY)

    # Create polygon
    poly_envelope = ogr.Geometry(ogr.wkbPolygon)
    poly_envelope.AddGeometry(ring)
    return poly_envelope

def getgrid(envelope):
    """
    Function with which to create the grid intersecting grid cells based
    on a feature geometry envelope. Initiate the class, and use getgrid() to
    generate a grid mesh and index information about the intersecting cells.

    getgrid() takes as input an OGR geometry envelope, and will
    compute the grid indices that contain the envelope.

    I-index values are numbered 0...n from the lower-left corner (left to right).
    J-index values are numbered 0...n from the lower-left corner (bottom to top).
    """

    # Calculate the number of grid cells necessary
    xmin, xmax, ymin, ymax = envelope

    # Find the i and j indices of the edges
    i0 = int((xmin - x_min) / cellsize // 1)  # Floor the value
    j0 = int((ymax - y_min) / cellsize // 1) + 1  # Floor the absolute value and add one
    i1 = int((xmax - x_min) / cellsize // 1) + 1  # Floor the value and add one
    j1 = int((ymin - y_min) / cellsize // 1)  # Floor the absolute value and add one

    # Calculate the rows, columns, and grid envelope
    cols = i1 - i0
    rows = j0 - j1
    gridEnvelope = ((i0 * cellsize) + x_min, (i1 * cellsize) + x_min, (j1 * cellsize) + y_min, (j0 * cellsize) + y_min)
    return cols, rows, gridEnvelope

def apply_calculation(shore_distance, max_depth=0, NoDataVal=-9999):
    """
    Apply distance calculation, which is each pixel's distance to shore, multiplied
    by the maximum depth, all divided by the maximum distance to shore. This provides
    a smooth slope from shore to lake max depth.

    shore_distance - Input numpy array containing array of distances to shoreline.
        Must contain positive values for distance away from shore and 0 elsewhere.
    max_depth - Input value with maximum lake depth.
    NoDataVal - value to assign to all non-lake pixels (zero distance to shore).
    """
    # Process the distance to shore and max depth into bathymetry
    bathymetry_arr = (shore_distance * max_depth) / float(shore_distance.max())
    bathymetry_arr[bathymetry_arr == 0] = NoDataVal
    return bathymetry_arr

# ---------- End Functions ---------- #

if __name__ == '__main__':
    tic = time.time()

    # Read the input maximum depth CSV
    df = pd.read_csv(lake_depths)
    depth = pd.Series(df[lake_depths_maxdepth_field].values, index=df[lake_depths_id_field]).to_dict()

    # Setup output and temporary file format using GDAL drivers
    OutrasterDriver = gdal.GetDriverByName(outFormat)  # GDAL raster driver for persistent rater files.
    if copyIntermediateToDisk:
        # In this case, the intermediate files will be written to disk rather than stay temporary
        TemprasterDriver = OutrasterDriver  # GDAL raster driver for persistent rater files.
    else:
        # In this case, the intermediate files will stay in memory.
        TemprasterDriver = gdal.GetDriverByName('Mem')  # GDAL memory driver for temporary raster files.
        dst_dataset1 = dst_dataset2 = ''  # Set dataset layer names to be blank
    TempVectorDriver = ogr.GetDriverByName('MEMORY')    # OGR vector driver for temporary vector files.

    # Open the lakes file in order to pull out the necessary basin geometries
    driver = ogr.GetDriverByName(inDriverName)  # Open the basin shapefile file with OGR, read-only access
    shp = driver.Open(inLakes, 0)  # 0 means read-only. 1 means writeable.

    # Make sure an object was returned when attempting to open the input vector file
    if shp is None:
        print("Open failed on {0}.\n".format(inLakes))
        raise SystemExit
    layer = shp.GetLayer()
    layerName = layer.GetName()
    basinProj = layer.GetSpatialRef()
    basinProjWKT = basinProj.ExportToWkt()
    numbasins = layer.GetFeatureCount()

    # Check for existence of provided fieldnames
    field_defn, fieldslist = checkfield(layer, fieldname, layerName)
    print('Number of input lakes: {0}'.format(numbasins))
    print('Projection: \n  {0}'.format(basinProj.ExportToPrettyWkt()))

    # If requested, initiate a shapefile where each feature will be the lake envelope
    if writeenvelopeSHP:
        # Create the shapefile that will store the subset raster polygons
        envDS = driver.CreateDataSource(envelopeSHP)
        envLayer = envDS.CreateLayer('', srs=basinProj, geom_type=ogr.wkbPolygon)
        if envLayer is None:
            print('        Envelope vector layer creation failed.\n')
            raise SystemExit
        envLayer.CreateField(field_defn)
        LayerDef = envLayer.GetLayerDefn()

    # Iterate over each geometry to calculate bathymetry
    tic2 = time.time()  # Counter for providing output every x iterations

    # To stop the iteration, provide a number for maxN < numbasins
    maxN = numbasins    # This will process all polygons in the input vector file
    i = 0
    for feat in layer:
        lk_id = feat.GetField(fieldname)

        # --- Prepare output layers --- #
        if copyIntermediateToDisk:
            dst_dataset1 = os.path.join(out_tiles_dir, '{0}_mask.{1}'.format(lk_id, suffix))
            dst_dataset2 = os.path.join(out_tiles_dir, '{0}_proximity.{1}'.format(lk_id, suffix))
        dst_dataset3 = os.path.join(out_tiles_dir, '{0}_bathymetry.{1}'.format(lk_id, suffix))

        # Get feature geometry and extent envelope
        geom = feat.GetGeometryRef()

        # Create temporary polygon vector layer from feature geometry so that we can rasterize it (Rasterize needs a layer)
        ds = TempVectorDriver.CreateDataSource('')
        Layer = ds.CreateLayer('', geom_type=ogr.wkbPolygon, srs=basinProj)  # Use projection from input vector layer
        outfeature = ogr.Feature(Layer.GetLayerDefn())      # Create the output feature
        outfeature.SetGeometry(geom)                        # Set geometry of output feature
        Layer.SetFeature(outfeature)

        # Get raster cell bounds (i,j) for polygon envelope bounds
        cols, rows, gridEnvelope = getgrid(geom.GetEnvelope())

        if writeenvelopeSHP:
            # The output rectangle polygon will represents each raster subset
            out_envelope = my_envelope(gridEnvelope)

            # Create feature for this rectangle and add it to the output layer
            feature = ogr.Feature(LayerDef)             # Create a new feature (attribute and geometry)
            feature.SetField(fieldname, lk_id)
            feature.SetGeometry(out_envelope.Clone())   # Make a feature from geometry object
            envLayer.CreateFeature(feature)
            feature.Destroy()
            feature = out_envelope = None

        # Build the first intermediate raster layer (mask array)
        mask_ds = TemprasterDriver.Create(dst_dataset1, cols, rows, 1, gdal.GDT_Byte)
        mask_ds.SetProjection(basinProjWKT)
        mask_ds.SetGeoTransform((gridEnvelope[0], cellsize, 0, gridEnvelope[3], 0, -cellsize))

        # Build additional raster layers and copy raster info to them (SpatialReference, Geotransform, etc)
        prox_ds = TemprasterDriver.Create(dst_dataset2, cols, rows, 1, gdal.GDT_Float32)

        # Add option for compression
        bathy_ds = OutrasterDriver.Create(dst_dataset3,
                                            cols,
                                            rows,
                                            1,
                                            gdal.GDT_Float32,
                                            options=['COMPRESS=LZW'])
        CopyDatasetInfo(mask_ds, prox_ds)
        CopyDatasetInfo(mask_ds, bathy_ds)

        # Rasterize the geometry
        gdal.RasterizeLayer(mask_ds, [1], Layer, None, None, [1], options=['ALL_TOUCHED=FALSE'])

        # Using target pixel values of 0. Distance units can be GEO or PIXEL. PIXEL is default.
        # proximity_options = ["VALUES=0", "DISTUNITS=GEO"]
        proximity_options = ["VALUES=0", "DISTUNITS=PIXEL"]
        srcband = mask_ds.GetRasterBand(1)
        dstband = prox_ds.GetRasterBand(1)
        stats = srcband.GetStatistics(0, 1)                                 # Calculate statistics
        gdal.ComputeProximity(srcband, dstband, options=proximity_options)
        stats = dstband.GetStatistics(0, 1)                                 # Calculate statistics

        # Process the distance to shore and max depth into bathymetry
        max_depth = depth[lk_id]
        bathymetry_arr = apply_calculation(dstband.ReadAsArray(), max_depth=max_depth, NoDataVal=NoDataVal)

        # Save the output and calculate statistics
        outband = bathy_ds.GetRasterBand(1)
        BandWriteArray(outband, bathymetry_arr)             # Write array to raster band
        outband.SetNoDataValue(NoDataVal)                   # Set NoData
        stats = outband.GetStatistics(0, 1)                 # Calculate statistics

        # Clean up intermediate GDAL/OGR objects by setting to None. Could also use .Destroy()
        Layer = outfeature = ds = mask_ds = prox_ds = srcband = dstband = bathy_ds = outband = geom = None
        del bathymetry_arr, stats, cols, rows, gridEnvelope, lk_id

        # Print messages to diagnose script speed
        i += 1
        if i % 1000 == 0:
            print('      1000 Records in {0:3.2f} seconds.'.format(time.time() - tic2))
            tic2 = time.time()
        # Stop after x number of features
        if i == maxN:
            break

    # Clean up and finish
    if envelopeSHP:
        envLayer = LayerDef = envDS = None
    feat = OutrasterDriver = TemprasterDriver = layer = None
    print('Created {0} raster tiles.'.format(i))
    del i, feat, geom
    print('Process completed in {0: 3.2f} seconds'.format((time.time() - tic)))