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


# Pixels to coordinates 
def pixelsToCoordinates(P1,P2):
    global transMatrix, transVector
    #Initial coordinates of 2 points (lower and lefter points: 1 and 2)
    p1 = {  'x': 129668.0,
            'y': -29854.0   }
    p2 = {  'x': 135757.0,
            'y': -14281.0   } 
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
    cos = V_norm[0][0]*v_norm[0][0]+V_norm[0][1]*v_norm[0][1]
    sin = math.sqrt(1-cos*cos)
    
    #Calculate the rotation matrix (P)* [(cos,-sin),(sin,cos)] = (p)
    res = np.array([[cos,-sin],[sin,cos]])
    transMatrix = np.linalg.inv(res)


# Pixels to coordinates 
def pixelsToCoordinates_old(P1,P2):
    global transMatrix, transVector
    #Initial coordinates of 2 points (lower and lefter points: 1 and 2)
#    p1 = {  'x': 430852,
#            'y': 4574491   }
#    p2 = {  'x': 432152,
#            'y': 4591231   }      
#    p1 = {  'x': 134385.0,
#            'y': -28817.0   }
    p1 = {  'x': 129668.0,
            'y': -29854.0   }
    p2 = {  'x': 135757.0,
            'y': -14281.0   } 
    #Calculate de translation vector, P1 - T = p1            
    Tx = P1[0] - p1['x']
    Ty = P1[1] - p1['y']
    transVector = np.array([Tx,Ty])
    
    #Transform the 2 provided points using the transVector
    P1 = P1 - transVector
    P2 = P2 - transVector
    
    #Calculate the transformation matrix (P)* [(a,b),(c,d)] = (p)
    a = (p1['y']*P2[0]-P1[0]*p2['y'])/(p1['y']*p2['x']-p1['x']*p2['y'])
    b = (p1['x']*P2[0]-P1[0]*p2['x'])/(p1['x']*p2['y']-p1['y']*p2['x'])
    c = (p1['y']*P2[1]-P1[1]*p2['y'])/(p1['y']*p2['x']-p1['x']*p2['y'])
    d = (p1['x']*P2[1]-P1[1]*p2['x'])/(p1['x']*p2['y']-p1['y']*p2['x'])
    #return the matrix values
    res = np.array([[a,b],[c,d]])
    transMatrix = np.linalg.inv(res)



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
            velArray = np.empty((width,height))
            for x in range(width):
                for y in range(height):
                    r, g, b = rgb_im.getpixel((x, y))
                    if (r, g, b) == (52, 52, 52):
                        if x > pix_max_x[0]:
                            pix_max_x = [x,y]
                        if x < pix_min_x[0]:
                            pix_min_x = [x,y]
                    rgb = (str(r), str(g), str(b))
                    velArray[x,y] = rgb_colors_values[', '.join(rgb)]
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
    # create the same mesh now over the raster geometry
    global meshCoordArray, meshCoordArrayResolution_x, meshCoordArrayResolution_y
    dim_X , dim_Y = geoArray.shape
    meshCoordArrayResolution_x = resolution_x * dim_x / cols
    meshCoordArrayResolution_y = resolution_y * dim_y / rows
    meshCoordArray = []
    for y in np.arange(meshCoordArrayResolution_y/2.0,resolution_y * dim_y, meshCoordArrayResolution_y):
        row_list = []
        for x in np.arange(meshCoordArrayResolution_x/2.0, resolution_x *dim_x, meshCoordArrayResolution_x):
            row_list.append((x, y))
        meshCoordArray.append(row_list)

    meshCoordArray = np.array(meshCoordArray)


    # create mesh over pixel 
    global meshPixelArray, meshPixelArrayResolution_x, meshPixelArrayResolution_y

    dim_X , dim_Y, dim_Z = meshCoordArray.shape
    meshPixelArray = np.zeros((dim_X , dim_Y))
    for Y in range(dim_Y):
        for X in range(dim_X):
            pixel = meshPixelArray[X][Y] - transVector
            meshCoordArray[X][Y] = transMatrix.dot(pixel)

    meshPixelArrayResolution_x = abs(meshPixelArray[0][0][0] -  meshPixelArray[0][1][0]) / 2.0
    meshPixelArrayResolution_y = abs(meshPixelArray[0][0][1] -  meshPixelArray[1][0][1]) / 2.0


""" Main Code """
readRasterData()
readInitialVelocityData(vel_from = 0, vel_to = 1)
pixelsToCoordinates(edgePoints[0]['bottom-left'], edgePoints[0]['top-rigth'])
createMesh(cols = 4, rows = 4)
np.absolute
