# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:37:39 2015
Script per calcular el recurs mini eolic:
Processat del mapa cartografic per generar un nou fitxer amb la info 
de cada punt del raster si pertany o no a un edifici, i si pertany a quin 
i quin tipus
@author: daniel
"""
# Libraries
from osgeo import gdal
import numpy as np
from numpy import nan
from numpy import isnan
import os.path
import math
import csv
import time
import copy
 


# Reading geometry data from the raster file (values are stored in a matrix called geoArray)
def readRasterData(file = 'barcelona_raster_augusto_500x400 v2.asc'):
    global geoArray, top_left_x, resolution_x, top_left_y, resolution_y
    #initial variables
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    ini_height = 386   # Value needed because the raster file is displaced in z
    if os.path.isfile(path_geo + file):
        print 'Reading initial geometry from file %s' % file
        src_ds = gdal.Open(path_geo + file)
        srcband = src_ds.GetRasterBand(1)  
        cols = src_ds.RasterXSize
        rows = src_ds.RasterYSize
        NoDataValue = srcband.GetNoDataValue()
        geoTransform = src_ds.GetGeoTransform()
        #print geoTransform
        top_left_x = geoTransform[0] # top left x 
        resolution_x = geoTransform[1] # w-e resolution 
        top_left_y = geoTransform[3] # top left y
        resolution_y = geoTransform[5] # n-s resolution 
        geoArray = srcband.ReadAsArray(0, 0, cols, rows)
        # Changing NoDataValue into dataArray to nan
        geoArray[geoArray == NoDataValue] = nan
        # Considering the initial_height
        #geoArray = geoArray - ini_height
    else:
        print 'file %s does not exist' % str(path_geo + file)


def modifyRasterFile(file = 'barcelona_raster_augusto_500x400 v5.asc'):
    #global rasterNewZValues
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    print('Calculating the new raster file')
    print('-------------------------------')
    if not os.path.isfile(path_geo + file):
        print('Calculating the coordinates of the raster file')
        dim_X , dim_Y = geoArray.shape
        rasterNewZValues = np.empty((dim_X,dim_Y)) 
        rasterNewZValues[:] = np.nan
        rasterCoordArray = []
        for y in np.arange(top_left_y, top_left_y+resolution_y*dim_Y, resolution_y):
            row_list = []
            for x in np.arange(top_left_x, top_left_x+resolution_x*dim_X, resolution_x):
                row_list.append((x, y))
            rasterCoordArray.append(row_list) 
        rasterCoordArray = np.array(rasterCoordArray)    
        print('Calculating the new position and new z')
        intervalPercentage = 5
        listPercentage = [int(dim_Y*x/100.0) for x in range(0,101,intervalPercentage)]
        t0 = time.time()
        listPoints_y = np.arange(top_left_y, top_left_y+resolution_y*dim_Y, resolution_y)  
        listPoints_y = listPoints_y.tolist()
        listPoints_x = np.arange(top_left_x, top_left_x+resolution_x*dim_X, resolution_x)
        listPoints_x = listPoints_x.tolist()
        for y in np.arange(dim_Y):
            for x in np.arange(dim_X):
                p = rasterCoordArray[y][x]
                new_p = modifyRasterPoint(p)
                #print p, new_p
                #new_p = [922098.7312161 , 4604435.92692865]
                # find the Square that contains the new point
                listPoints_y0 = copy.copy(listPoints_y)
                listPoints_y0.append(new_p[1])
                listPoints_y0.sort()
                new_p_y_position = listPoints_y0.index(new_p[1])
                y0 = listPoints_y0[new_p_y_position-1] if new_p_y_position > 0 else nan
                y1 = listPoints_y0[new_p_y_position+1] if new_p_y_position < dim_Y else nan
                listPoints_x0 = copy.copy(listPoints_x)
                listPoints_x0.append(new_p[0])
                listPoints_x0.sort()
                new_p_x_position = listPoints_x0.index(new_p[0])
                x0 = listPoints_x0[new_p_x_position-1] if new_p_x_position > 0 else nan
                x1 = listPoints_x0[new_p_x_position+1] if new_p_x_position < dim_X else nan
                square = [[x0,y0],[x0,y1],[x1,y1],[x1,y0]]
                #print square
                #Calculating the weigth of the points of the square
                dist = []
                for point in square:
                    v = np.array(new_p) - np.array(point)
                    d = np.sqrt(v.dot(v))
                    inv_dist = (1 / d) if d != 0 else 1e6
                    dist.append(inv_dist)
                dist = np.array(dist)
                mdist = np.ma.masked_array(dist,np.isnan(dist))
                mdistSum = np.sum(mdist)
                weigth = dist / mdistSum
                #print weigth
                #Calculating the new z
                list_z = []
                for point in square:
                    if point[1] is not nan and point[0] is not nan:
                        p_y_position = listPoints_y.index(point[1]) 
                        p_x_position = listPoints_x.index(point[0])
                        if p_y_position < dim_Y and p_x_position < dim_X :
                            #print p_x_position, p_y_position
                            new_z = geoArray[p_y_position][p_x_position]
                        else:
                            new_z = np.nan
                    else:
                        new_z = np.nan
                    list_z.append(new_z)
                #print list_z
                list_z = np.array(list_z)
                z_weigth = list_z * weigth
                #print z_weigth
                mZ = np.ma.masked_array(z_weigth,np.isnan(z_weigth))
                #print mZ
                mZSum = np.sum(mZ)
                #print mZSum
                rasterNewZValues[y][x] = mZSum
                break
            if y in listPercentage:
                print '%s%% done in %s seconds' % (str(y*100/dim_Y), str(round(time.time()-t0,0)))
                t0 = time.time()

        print('Saving the new raster file')
        fl = open(path_geo + file, 'w')
        fl.write('NCOLS %s \n' %dim_X)
        fl.write('NROWS %s \n' %dim_Y)
        fl.write('XLLCORNER %s \n' %top_left_x)
        fl.write('YLLCORNER %s \n' %(top_left_y+resolution_y*dim_Y))
        fl.write('CELLSIZE %s \n' %resolution_x)
        fl.write('NODATA_VALUE -9999 \n')
        rasterNewZValues[np.isnan(rasterNewZValues)] = int(-9999)
        writer = csv.writer(fl, delimiter=' ')
        for values in rasterNewZValues:
            writer.writerow(values)
        fl.close()    
    else:    
        print 'file %s already exists' % str(path_geo + file)



# Coordinates to pixels
def modifyRasterPoint(p):
    #Ini variables
    # This ini variables were used to make a first modification
    #ang = 3.8*math.pi/180.0 
    #origin_xy = [933208, 4603698]  
    #This second is needed to full fill the transformation
    ang = 0.2*math.pi/180.0 
    origin_xy = [931624, 4596432] 
    #scale = 1
    #translationVector=[0,0]
    
    #Change the point to an array
    if not isinstance(p, np.ndarray):
        p = np.array(p)
    
    #Find the vector p - origin    
    v = np.array(p) - np.array(origin_xy)
    
    #Calculate the module
#    mod_v = np.sqrt(v.dot(v))
    
    #Applying the rotation matrix
    transfMatrix = np.array([[math.cos(ang),math.sin(ang)],[-math.sin(ang),math.cos(ang)]])
    v_rotated =  transfMatrix.dot(v)
    
    #Applying the scale factor
#    if mod_v != 0 :
#        mod_scaled_v = v_rotated * scale
#    else:
#        mod_scaled_v = v_rotated
    
    #The scaled point
#    mod_scaled_p = mod_scaled_v + origin_xy
    
    #Applying the translation            
#    p_new = mod_scaled_p + translationVector

    p_new = v_rotated + origin_xy
    
    #Return the new value
    return p_new    




#Return the distance between two points
def calculateDistance(p, centralPoint):
    if not isinstance(p, np.ndarray):
        p = np.array(p)
    if not isinstance(centralPoint, np.ndarray):
        centralPoint = np.array(centralPoint)
    v = p - centralPoint
    dist = np.sqrt(v.dot(v))
    return dist
    
    


        
    

""" Main Code """
readRasterData(file = 'barcelona_raster_augusto_3000x3000 v8.asc')
modifyRasterFile(file = 'barcelona_raster_augusto_3000x3000 v10.asc')


