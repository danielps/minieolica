# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:37:39 2015

Script per calcular el recurs mini eolic

@author: daniel
"""
# Libraries
from PIL import Image
from osgeo import gdal
import numpy as np
import os.path
import math


# Find the different information (translation, rotation and scalefactor) needed to change pixels to coordinates 
def pixelsToCoordinatesIni(P1,P2):
    global transMatrix, transVector, scaleFactor, originPointPixel, originPointCoord
    #Initial coordinates of 2 points (lower and lefter points: 1 and 2)
    print 'Calculating the needed information (translation, rotation and scalefactor) to change coordinates-raster to pixel-velocity'
    p1 = {  'x': 129515.0,
            'y': -29683.0   }
    p2 = {  'x': 135770.0,
            'y': -14263.0   } 
    #Calculate de translation vector, P1 - T = p1            
    Tx = P1[0] - p1['x']
    Ty = P1[1] - p1['y']
    transVector = np.array([Tx,Ty])
       
    #Calculate the 2 vectors
    V = np.array(P2) - np.array(P1)
    v = np.array([p2['x']-p1['x'], p2['y']-p1['y']])
    
    #Calculate the 2 modules
    mod_V = np.sqrt(V.dot(V))
    mod_v = np.sqrt(v.dot(v))
    
    #Calculate the 2 norm vectors
    V_norm = V / mod_V
    v_norm = v / mod_v

    #Calculate the angle between them
    coseno = V_norm[0]*v_norm[0]+V_norm[1]*v_norm[1]
    seno = math.sqrt(1-coseno*coseno)
    
    #Calculate the rotation matrix (P)* [(cos,-sin),(sin,cos)] = (p)
    res = np.array([[coseno,-seno],[seno,coseno]])
    transMatrix = np.linalg.inv(res)
    
    #Calculate the scale factor, and the origin point to apply the scaling
    scaleFactor = mod_V / mod_v
    originPointPixel = np.array(P1)
    originPointCoord = np.array([p1['x'],p1['y']])
    

# Coordinates to pixels
def coordinatesToPixels(p):
#    p = [124660,-18945]
#    p = [129668,-29854]
#    p = [135757,-14281]
    #Change the point to an array
    if not isinstance(p, ndarray):
        p = np.array(p)
    
    #Find the vector p - origin    
    v = np.array(p) - np.array(originPointCoord)
    
    #Calculate the module
    mod_v = np.sqrt(v.dot(v))
    
    #Applying the rotation matrix
    #transMatrix_0 = np.linalg.inv(transMatrix)
    v_rotated =  transMatrix.dot(v)
    
    #Applying the scale factor
    if mod_v != 0 :
        mod_scaled_v = v_rotated * scaleFactor
    else:
        mod_scaled_v = v_rotated
    
    #The scaled point
    mod_scaled_p = mod_scaled_v + originPointCoord
    
    #Applying the translation            
    p_new = mod_scaled_p + transVector

    #Return the new value
    return p_new
    

# Reading velocity data (values are stored in a matrix called velArray)
def readInitialVelocityData(vel_from = 0, vel_to = 39):
    global velArray, edgePoints
    #rgb to number of yearly hours
    initial = 250
    step = 250
    rgb_colors_values = {
        '188, 249, 253': 0,
        '163, 225, 249': initial + 1*step,
        '135, 201, 247': initial + 2*step,
        '107, 183, 244': initial + 3*step,
        '78, 163, 241' : initial + 4*step,
        '44, 141, 238' : initial + 5*step,
        '23, 117, 227' : initial + 6*step,
        '23, 91, 213'  : initial + 7*step,
        '21, 70, 196'  : initial + 8*step,
        '17, 46, 178'  : initial + 9*step,
        '9, 23, 165'   : initial + 10*step,
        '0, 0, 147'    : initial + 11*step,
        '52, 52, 52'   : None,      #black
        '255, 255, 255': None       #white
    }
    edgePoints = {}
    path_vel = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Dades inicials/Mapes minieolica/'
    for i in range(vel_from,vel_to):
        seq = (str(i),str(i+1))
        file = '-'.join(seq) + '.bmp'
        if os.path.isfile(path_vel + file):
            print 'Reading initial velocity from file %s' % file
            im = Image.open(path_vel + file)
            (width, height) = im.size
            width = int(width)
            width_25 = int(width*0.25)
            width_50 = int(width*0.5)
            width_75 = int(width*0.75)
            width_100 = int(width-1)
            height = int(height)
            rgb_im = im.convert('RGB')       
            # Insert the values into a 2D matrix (x= width, y=height)
            # and look for the reference points need to pixelToCoordinates
            pix_max_x = [0,0]
            pix_min_x = [width, height]
            velArray = np.empty((width,height+1))
            for x in range(width):
                for y in range(height):
                    r, g, b = rgb_im.getpixel((x, y))
                    if (r, g, b) == (52, 52, 52):
                        if x > pix_max_x[0]:
                            pix_max_x = [x,abs(y-height)]
                        if x < pix_min_x[0]:
                            pix_min_x = [x,abs(y-height)]
                    rgb = (str(r), str(g), str(b))
                    velArray[x,abs(y-height)] = rgb_colors_values[', '.join(rgb)]
                if (x == width_25) : print '25% done'
                if (x == width_50) : print '50% done'
                if (x == width_75) : print '75% done'
                if (x == width_100): print '100% done'
            edgePoints[i] = {'top-rigth': pix_max_x, 'bottom-left': pix_min_x}
            

# Reading geometry data from the raster file (values are stored in a matrix called geoArray)
def readRasterData(file = 'barcelona_raster_augusto_500x400.asc'):
    global geoArray, top_left_x, resolution_x, top_left_y, resolution_y
    #rgb to number of yearly hours
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    if os.path.isfile(path_geo + file):
        print 'Reading initial geometry from file %s' % file
        src_ds = gdal.Open(path_geo + file)
        srcband = src_ds.GetRasterBand(1)  
        cols = src_ds.RasterXSize
        rows = src_ds.RasterYSize
        NoDataValue = srcband.GetNoDataValue()
        geoTransform = src_ds.GetGeoTransform()
        top_left_x = geoTransform[0] # top left x 
        resolution_x = geoTransform[1] # w-e pixel resolution 
        top_left_y = geoTransform[3] # top left y
        resolution_y = geoTransform[5] # n-s pixel resolution 
        geoArray = srcband.ReadAsArray(0, 0, cols, rows)
        # Changing NoDataValue into dataArray to nan
        geoArray[geoArray == NoDataValue] = nan
        #How to calculate the mean excluding nan values
        mdat = np.ma.masked_array(geoArray[1,],np.isnan(geoArray[1,]))
        mm = np.mean(mdat)
        mm1 = np.sum(mdat)
        mm2 = np.std(mdat)
        mm3 = np.min(mdat)
        mm4 = np.max(mdat)
    else:
        print 'file %s does not exist' % str(path_geo + file)

# Create the mesh of point where the calculations will be donemeshCoordArray
# Resolution is 1/2 of the distance between two adjacent points
def createMesh(cols = 400, rows = 400):
    # create a mesh over the raster geometry
    global meshCoordArray, meshCoordArrayResolution_x, meshCoordArrayResolution_y
    print 'Creating mesh over Raster geometry'
    dim_X , dim_Y = geoArray.shape
    meshCoordArrayResolution_x = resolution_x * dim_X / cols
    meshCoordArrayResolution_y = resolution_y * dim_Y / rows
    meshCoordArray = []
    for y in np.arange(top_left_y+meshCoordArrayResolution_y/2.0,top_left_y+resolution_y*dim_Y, meshCoordArrayResolution_y):
        row_list = []
        for x in np.arange(top_left_x+meshCoordArrayResolution_x/2.0, top_left_x+resolution_x*dim_X, meshCoordArrayResolution_x):
            row_list.append((x, y))
        meshCoordArray.append(row_list) 

    meshCoordArray = np.array(meshCoordArray)
    
    # create the same mesh now over the over pixel file
    print 'Creating mesh over velocity pixel file'
    global meshPixelArray, meshPixelArrayResolution_x, meshPixelArrayResolution_y

    dim_X , dim_Y, dim_Z = meshCoordArray.shape
    meshPixelArray = []
    for Y in range(dim_Y):
        row_list = []
        for X in range(dim_X):
            #We only need to change the coordinates to pixels and store the new value in the same position
            pixel = coordinatesToPixels(meshCoordArray[X][Y])
            row_list.append(pixel.tolist())
        meshPixelArray.append(row_list)

    meshPixelArray = np.array(meshPixelArray)

    meshPixelArrayResolution_x = abs(meshPixelArray[0][0][0] -  meshPixelArray[0][1][0]) / 2.0
    meshPixelArrayResolution_y = abs(meshPixelArray[0][0][1] -  meshPixelArray[1][0][1]) / 2.0

#Calculate the mean velocity value over the pixel mesh
def calculateMeanVel():
    global meshMeanVelArray
    print('Calculating the mean velocities')
    dim_X , dim_Y, dim_Z = meshPixelArray.shape
    meshMeanVelArray = []
    for y in np.arange(dim_Y):
        row_list = []
        for x in np.arange(dim_X):
            # Find the edges of the square in Coord and
            central_point_coord = meshCoordArray[x][y]
            point_top_left_coord = central_point_coord + np.array([-abs(meshCoordArrayResolution_x), abs(meshCoordArrayResolution_y)])
            point_bot_left_coord = central_point_coord + np.array([-abs(meshCoordArrayResolution_x), -abs(meshCoordArrayResolution_y)])
            point_top_righ_coord = central_point_coord + np.array([abs(meshCoordArrayResolution_x), abs(meshCoordArrayResolution_y)])
            point_bot_righ_coord = central_point_coord + np.array([abs(meshCoordArrayResolution_x), -abs(meshCoordArrayResolution_y)])
            # then transform it into Pixel
            central_point_pixel = meshPixelArray[x][y]
            point_top_left_pixel = coordinatesToPixels(point_top_left_coord)
            point_bot_left_pixel = coordinatesToPixels(point_bot_left_coord)
            point_top_righ_pixel = coordinatesToPixels(point_top_righ_coord)
            point_bot_righ_pixel = coordinatesToPixels(point_bot_righ_coord)
            # Look for the points inside this square
            
            row_list.append((x, y))
        meshCoordArray.append(row_list) 

    meshCoordArray = np.array(meshCoordArray)
    

""" Main Code """
readRasterData()
readInitialVelocityData(vel_from = 0, vel_to = 1)
pixelsToCoordinatesIni(edgePoints[0]['bottom-left'], edgePoints[0]['top-rigth'])
createMesh(cols = 4, rows = 4)


