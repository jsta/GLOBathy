# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Name:         WGS_84_cell_dimesion_calculator
# Purpose:      This script calculates cell dimensions based on latitude in WGS84 coordinate system.
#               There are two options:
#                   1) latitude for a single point can be entered in script inputs (below)
#                   2) program reads geocoded raster files from a path defined in script inputs (below)
#               Results of option (1) are printed to screen and option (2) are written in excel files.
#               Important Note:
#                   * Valid range for latitude is [-90, 90] degrees
# Created:      2020
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

# import modules 
globals().clear()
import os
import math
import glob
import time
import numpy as np
import pandas as pd
from osgeo import gdal,ogr

# Start timer
tic = time.time()

############################## [ Script Inputs] ###############################
# Please complete this section before running the script:

# Define if your raster input is not geocoded or you need 
# calculations for a single point. Valid range is [-90.0, 90.0]
latitude = 45.00

# Define path to your bathymetry files, otherwise leave empty
# and enter average latitude above
bathy_path = '/path/to/file(s).tif'

# Image resolution in arc-seconds (For GLOBathy use 1.0)
resolution = 1.0    
########################### [ End of Script Inputs] ###########################

############################# [ Script Functions] #############################
# define a function of that estimates cell dimensions based on latitude

#2020
#WGS84 equatorial radius (a): 6,378,137.00
#WGS84 polar radius (b): 6,356,752.30
#perimeter of WGS84 ellipse: 40,007,862.87

a = 6378137.00                      # WGS84 Equatorial Radius
b = 6356752.30                      # WGS84 Polar Radius
degrees = resolution/360.0/3600.0   # Distance to measure in fractions of a perimeter (One arc-second)

# Perimeter of the circle with radius 6378137 (equator - WGS84)
Equator_Circumference = (a * 2.0) * math.pi
print('Circumference of the Equator (m): {0}'.format(Equator_Circumference))

# https://www.fws.gov/r7/nwr/Realty/data/LandMappers/Public/Help/HTML/R7-Public-Land-Mapper-Help.html?Degreesandgrounddistance.html
Elipse_Perimeter = math.pi * ((3.0*(a+b))-math.sqrt((a+(3.0*b))*((3.0*a)+b)))
print('Perimeter of the elipse (N-S direction): {0} m'.format(Elipse_Perimeter))
    
def func_dim(latitude):
    lat = latitude                      # Latitude of measurement

    # Based on WGS84, the datum we are using if I understand correctly, I show c1_lat = 30.818404825 c1_lon =  30.922
    c1 = Elipse_Perimeter * degrees                                                 # North-South direction distance
    c2 = (Equator_Circumference * degrees) * math.cos(lat*math.pi/180.0)            # East-West direction distance at specified latitude
    return c1, c2
######################### [ End of Script Functions] ##########################

# Start calculation of cell dimensions:
if len(bathy_path)==0:      # Single point calculation
    if latitude<=90.0 and latitude>=-90.0:
        c_NS,c_EW = func_dim(latitude)
        print('North-South distance of each cell at any latitude: {0} m'.format(c_NS))
        print('East-West distance of each cell at latitude {0}: {1} m'.format(latitude, c_EW))
    else:
        print('Error: Enter a latitude between -90.0 and 90.0!')
        print('Program is terminated !!!')
            
else:                       # Raster calculation
    bathyfiles=glob.glob(bathy_path+'*.tif')
    NUMfiles=len(bathyfiles)
    
    for i in range(0,NUMfiles):
        string = "{:8.0f}: Now Processing ---> File {:20s} ...... {:6.2f}% completed in {:8.1f} seconds !".format(i+1, str(bathyfiles[i].replace(bathy_path, '')), ((i+1) / NUMfiles * 100), time.time()-tic)
        print(string)
        
        # read geocoded raster
        ds = gdal.Open(bathyfiles[i])
        width = ds.RasterXSize
        height = ds.RasterYSize
        gt = ds.GetGeoTransform()
        minx = gt[0]      # min longitude
        maxx = gt[0] + width*gt[1] + height*gt[2]     # max longitude
        miny = gt[3] + width*gt[4] + height*gt[5]     # min latitude
        maxy = gt[3]        # max latitude
        lat_center = (miny+maxy)/2.0      # average latitude
        
        if maxy<=90.0 and maxy >=-90.0 and miny<=90.0 and miny >=-90.0:
            c_temp,c_EW_avg = func_dim(lat_center) # average EW dimension
            print('          North-South distance of each cell at any latitude: {0} m'.format(c_temp))
            print('          Average East-West distance of each cell at average latitude {0}: {1} m'.format(lat_center, c_EW_avg))
            
            # generate the latitude 2D array
            c0 = np.linspace(maxy, miny, num=height)
            latitude = np.transpose([c0] * width)

            # estimate 2D cell dimension
            c_NS, c_EW = np.vectorize(func_dim)(latitude)

            # dataframe for converted cell dimension:
            df1 = pd.DataFrame(c_EW)

            # dataframe for converted other dimensions:
            other_info = {'Direction': ['N-S (the entire domian)','E-W (avarge in the entire domain)','','',
                                        'Min Latitude','Max Latitude','Min Longitude','Max Longitude'],
                          'Cell Dimension (m)': [c_NS[0,0], c_EW_avg, '','',
                                                 miny,maxy,minx,maxx]}
            df2 = pd.DataFrame(other_info, columns = ['Direction', 'Cell Dimension (m)'])
        
            # write data to excel:
            with pd.ExcelWriter(str(bathyfiles[i].replace('.tif', '.xlsx'))) as writer:
                df1.to_excel(writer, sheet_name='East_West 2D Cell Dimension (m)')
                df2.to_excel(writer, sheet_name='Other Dimensions and Info', index=False)
            
        else:
             print('          Error: Only latitude in range -90.0 and 90.0 is valid!')
             print('          Current file is skipped !!!')