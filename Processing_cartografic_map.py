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
from osgeo import gdal_array
from osgeo import osr
import numpy as np
from numpy import nan
import os.path
import math
import csv
from osgeo.gdalconst import *
import time

 

# Reading the cartographic map. The info is stored into cartographicMapArray(parcelId,RiskueId,ListCoorXY,center)
def readInitialMap(filename='Predios_bcn_ETRS89.srl'):
    global cartographicMapArray
    t0 = time.time()
    cartographicMapArray = []
    path_cad = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Cadastro/'
    #path_cad = '/home/xcipriano/Escriptori/MiniEolica/Cadastro/'
    if os.path.isfile(path_cad + filename):
        print 'Reading serialized cadastre from file %s' % filename
        file = open(path_cad + filename, 'r')
        #i = 0
        for line in file:
            #i += 1
            #Seeking the info related to coordinates that define the polygon
            if line.find(';') >= 0:
                coordinates = line.split(';')[-1]
                coordinates = coordinates.split('|')[0:-1]
                #print 'coor: ', coordinates
                newCoor = [[float(coordinates[n]), float(coordinates[n+1])] for n in xrange(0,len(coordinates),2)]
                center = np.mean(newCoor, axis=0)
                min_xy = np.min(newCoor, axis=0)
                max_xy = np.max(newCoor, axis=0)
                square = [[min_xy[0], max_xy[1]], [max_xy[0], max_xy[1]], [max_xy[0], min_xy[1]], [min_xy[0], min_xy[1]], [min_xy[0], max_xy[1]]]
                v = max_xy - min_xy
                dist = np.sqrt(v.dot(v))
                #print newCoor
                #print 'min: ', min_xy
                #print 'max: ', max_xy
                #print 'square: ', square
                #print 'center: ', center
                #print 'dist: ', dist
            else:
                parcelId = line.split('|')[4]
                parcelId = parcelId.replace('#RISKUE', '')
                RiskueId = line.rstrip().split('|')[-1]
                #print 'parcel: ', parcelId
                #print 'riskue: ', RiskueId
                if RiskueId != '0' :
                    newCoor = readInitialMapCheckPoints(newCoor)
                    cartographicMapArray.append([parcelId,RiskueId,newCoor,center.tolist(),square,dist/2.0])
            #if i > 3: break
    else: 
        print 'file %s does not exist' % str(path_cad + filename)      
    
    print 'Reading serialized cadastre has finished in %s seconds' % str(round(time.time()-t0,1))

#Checking the list of point and return only the needed points (to avoid points closer to another or points over the same line)
def readInitialMapCheckPoints(listPoints):
    import copy
    # Eliminate the points over the same line
    newListPoints = copy.deepcopy(listPoints)
    for i in range(0,len(listPoints)-2):     
        area = triangleArea(listPoints[i], listPoints[i+1], listPoints[i+2])
        #print listPoints[i], listPoints[i+1], listPoints[i+2], area
        if abs(area) < 1e-1:
            #print 'remove 2 ', listPoints[i+1]
            newListPoints.remove(listPoints[i+1])
    #print newListPoints
    #Eliminate the points that are closer
    newListPoints2 = copy.deepcopy(newListPoints)
    for i in newListPoints2:
        for j in newListPoints2:
            #print i, ' ', j#, ' ', calculateDistance(i,j)
            if i != j and calculateDistance(i,j) < 3.0:
                #print 'removed point: %s' % j
                newListPoints2.remove(j)
    # To ensure that the polygon is closed
    if newListPoints2[0] != newListPoints2[-1]:
        newListPoints2.append(newListPoints2[0])
    #print newListPoints2
    if len(newListPoints2) >=5 :
        return newListPoints2   
    else:
        return listPoints


#Creating a region where the different parcels will be included
def readInitialMapRegionClassified(n_elements = 10):
    global cartographicMapRegion, mapRegion, cartographicMapRegionIndex

    dim_X , dim_Y = geoArray.shape
    # np.arange does not include the stop (2on term), because of this I have to increase the stop 
    list_y = np.arange(top_left_y, -1 + top_left_y+resolution_y*dim_Y, resolution_y*dim_Y/n_elements)
    list_y = list_y.tolist()
    list_x = np.arange(top_left_x, 1 + top_left_x+resolution_x*dim_X, resolution_x*dim_X/n_elements)
    list_x = list_x.tolist()
    cartographicMapRegion = np.zeros((n_elements,n_elements), dtype=object)     
    cartographicMapRegionIndex = np.zeros((n_elements,n_elements), dtype=object)            
    mapRegion = np.zeros((n_elements,n_elements), dtype=object)
    for y in range(1,len(list_y)):
        y_0 = list_y[y-1]
        y_1 = list_y[y]
        for x in range(1,len(list_x)):
            print y-1, x-1
            x_0 = list_x[x-1]
            x_1 = list_x[x]
            parcel_inside = []
            parcel_index = []
            square = [[x_0,y_0],[x_0,y_1],[x_1,y_1],[x_1,y_0],[x_0,y_0]]
            d = np.array([[x_0,y_0],[x_0,y_1],[x_1,y_1],[x_1,y_0]])
            squareCenter = np.mean(d, axis=0)
            min_xy = np.min(d, axis=0)
            max_xy = np.max(d, axis=0)
            v = max_xy - min_xy
            squareDist = np.sqrt(v.dot(v))
            #print square
            i = 0
            for parcel in cartographicMapArray:
                square_to_check = parcel[4]
                parcelDist = parcel[5]
                parcelCenter = parcel[3]
                dist = calculateDistance(squareCenter, parcelCenter)
                if dist < parcelDist + squareDist:
                    for point in square_to_check:
                        isInsideSquare = WindingNumber(point, square)
                        if abs(isInsideSquare) > 5e-2:
                            parcel_inside.append(parcel[0])
                            parcel_index.append(i)
                            break
                i += 1
            cartographicMapRegion[x-1][y-1] = parcel_inside
            cartographicMapRegionIndex[x-1][y-1] = parcel_index
            mapRegion[x-1][y-1] = square
        
        
# Reading geometry data from the raster file (values are stored in a matrix called geoArray)
def readRasterData(file = 'barcelona_raster_augusto_500x400 v2.asc'):
    global geoArray, top_left_x, resolution_x, top_left_y, resolution_y
    #initial variables
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    #ini_height = 386   # Value needed because the raster file is displaced in z
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



# Calculating the solid angle of a point respect to one surface. If the solid angle is 0, that means that the
# point is outside, if the value is 2PI the point is inside, and any other value the point is over one of the 
# edges of the surface.
def WindingNumber (p, ListSurfacePoints):
    totalAngle = 0
    qx = p[0]
    qy = p[1] 
    for i in range(1,len(ListSurfacePoints)):
        p1 = ListSurfacePoints[i-1]
        p2 = ListSurfacePoints[i]
        cx_p1 = p1[0]
        cy_p1 = p1[1]
        cx_p2 = p2[0]
        cy_p2 = p2[1] 
        # vectors
        u1 = qx - cx_p1
        v1 = qy - cy_p1
        u2 = qx - cx_p2
        v2 = qy - cy_p2    	
        # Calculating the producto escalar
        mod = math.sqrt(pow(u1,2)+pow(v1,2)) * math.sqrt(pow(u2,2)+pow(v2,2))
        if mod != 0:
            prod_escalar = (u1*u2 + v1*v2) / mod
            # To avoid rounding errors  
            if prod_escalar > 1:
                prod_escalar = 1
            if prod_escalar < -1:
                prod_escalar = -1            
            # Calculo el signo del producto vectorial a través del area del triangulo que forman los 3 puntos
            prod_vect = qx*cy_p1 - qy*cx_p1 + qy*cx_p2 - qx*cy_p2 + cx_p1*cy_p2 - cy_p1*cx_p2
            # Calculo el ángulo
            ang_sen = 0
            ang_cos = math.acos(prod_escalar)
            if prod_vect > 1e-6:
                ang_sen = 1 
            if prod_vect < -1e-6: 
                ang_sen = -1
                
            #print '%s %s %s %s %s %s'  % (p1, p2, mod, prod_escalar, prod_vect, totalAngle)
            totalAngle = totalAngle + ang_cos*ang_sen
         
    return totalAngle


#Calculate the triangle area
def triangleArea (p, p1, p2):
    #t0 = time.time()
    qx = p[0]
    qy = p[1] 
    cx_p1 = p1[0]
    cy_p1 = p1[1]
    cx_p2 = p2[0]
    cy_p2 = p2[1] 
    # Calculo el area del triangulo que forman los 3 puntos
    area = qx*cy_p1 - qy*cx_p1 + qy*cx_p2 - qx*cy_p2 + cx_p1*cy_p2 - cy_p1*cx_p2  
    #print time.time()-t0
    return area


#Return the distance between two points
def calculateDistance(p, centralPoint):
    if not isinstance(p, np.ndarray):
        p = np.array(p)
    if not isinstance(centralPoint, np.ndarray):
        centralPoint = np.array(centralPoint)
    v = p - centralPoint
    dist = np.sqrt(v.dot(v))
    return dist
    
    

#Calculate the mean velocity value over the pixel mesh
def calculateRasterInsideCartographic():
    global rasterInsideCartographicArray
    print('Calculating the coordinates of the raster file')
    dim_X , dim_Y = geoArray.shape
    rasterInsideCartographicArray = np.zeros((dim_X,dim_Y), dtype=object) 
    rasterCoordArray = []
    for y in np.arange(top_left_y, top_left_y+resolution_y*dim_Y, resolution_y):
        row_list = []
        for x in np.arange(top_left_x, top_left_x+resolution_x*dim_X, resolution_x):
            row_list.append((x, y))
        rasterCoordArray.append(row_list) 

    rasterCoordArray = np.array(rasterCoordArray)    

    print('Calculating the raster points inside cartographic')
    width_01 = int(dim_Y*0.01)
    width_05 = int(dim_Y*0.05)
    width_10 = int(dim_Y*0.1)
    width_20 = int(dim_Y*0.2)
    width_30 = int(dim_Y*0.3)
    width_40 = int(dim_Y*0.4) 
    width_50 = int(dim_Y*0.5)
    width_60 = int(dim_Y*0.6)
    width_70 = int(dim_Y*0.7)
    width_80 = int(dim_Y*0.8)
    width_90 = int(dim_Y*0.9)
    width_100 = int(dim_Y-1)
    minDistance = 1500
    t0 = time.time()
    for y in np.arange(dim_Y):
        for x in np.arange(dim_X):
            #print x, ' ', y
            # Find the edges of the square in Coord and
            p = rasterCoordArray[y][x]
            #print p
            #then transform it into UTM coordinates (not needed because the raster file is now georeferenced)
            #p = rasterToCoordinates(p)   
            #print p
            #p = [932952.591045, 4601664.363319]
            #p = [932814, 4601531]
            #and then I look for the inside points and points over the edges using the windingNumber
            for i in range(0,len(cartographicMapArray)):
                center = cartographicMapArray[i][3]
                distance = calculateDistance(center, p)
                if distance <= minDistance:
                    isInsideSquare = WindingNumber(p, cartographicMapArray[i][4])
                    if abs(isInsideSquare) > 5e-2:
                        #print cartographicMapArray[i][0]
                        listSurfacePoints = cartographicMapArray[i][2]
                        isInside = WindingNumber(p, listSurfacePoints)
                        if abs(isInside) > 5e-1:
                            #print 'is inside ', isInside
                            ini = rasterInsideCartographicArray[y][x]
                            if ini == 0:
                                rasterInsideCartographicArray[y][x] = [[cartographicMapArray[i][0],cartographicMapArray[i][1]]]
                            else: 
                                ini.append([cartographicMapArray[i][0],cartographicMapArray[i][1]])
                                rasterInsideCartographicArray[y][x] = ini
            #print rasterInsideCartographicArray[y][x]            
        if (y == width_01) : 
            print '1%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
            break
        if (y == width_05) : 
            print '5%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_10) : 
            print '10%% donein %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_20) : 
            print '20%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_30) : 
            print '30%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_40) : 
            print '40%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_50) : 
            print '50%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_60) :
            print '60%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_70) : 
            print '70%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_80) : 
            print '80%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_90) : 
            print '90%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_100): 
            print '100%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        break
    print '1 line done in %s seconds' % str(round(time.time()-t0,0))
    print('Saving the info cartographic into file')
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    fl = open(path_geo + "infoCartographic.csv", 'w')
    writer = csv.writer(fl, delimiter=';')
    for values in rasterInsideCartographicArray:
        writer.writerow(values)
    fl.close()    

   

#Calculate the mean velocity value over the pixel mesh
def calculateRasterInsideCartographic2():
    global rasterInsideCartographicArray
    print('Calculating the coordinates of the raster file')
    dim_X , dim_Y = geoArray.shape
    rasterInsideCartographicArray = np.zeros((dim_X,dim_Y), dtype=object) 
    rasterCoordArray = []
    for y in np.arange(top_left_y, top_left_y+resolution_y*dim_Y, resolution_y):
        row_list = []
        for x in np.arange(top_left_x, top_left_x+resolution_x*dim_X, resolution_x):
            row_list.append((x, y))
        rasterCoordArray.append(row_list) 

    rasterCoordArray = np.array(rasterCoordArray)    
    print('Calculating the raster points inside cartographic')
    width_01 = int(dim_Y*0.01)
    width_05 = int(dim_Y*0.05)
    width_10 = int(dim_Y*0.1)
    width_20 = int(dim_Y*0.2)
    width_30 = int(dim_Y*0.3)
    width_40 = int(dim_Y*0.4) 
    width_50 = int(dim_Y*0.5)
    width_60 = int(dim_Y*0.6)
    width_70 = int(dim_Y*0.7)
    width_80 = int(dim_Y*0.8)
    width_90 = int(dim_Y*0.9)
    width_100 = int(dim_Y-1)
    #minDistance = 1500
    t0 = time.time()
    for y in np.arange(dim_Y):
        for x in np.arange(dim_X):
            #print x, ' ', y
            # Find the edges of the square in Coord and
            p = rasterCoordArray[y][x]
            #print p
            #and then I look for the inside points and points over the edges using the windingNumber
            for i in range(0,len(cartographicMapArray)):
                center = cartographicMapArray[i][3]
                distance = calculateDistance(center, p)
                minDistance = cartographicMapArray[i][5]
                if distance <= minDistance:
                    #isInsideSquare = WindingNumber(p, cartographicMapArray[i][4])
                    isInsideSquare = 1
                    if abs(isInsideSquare) > 5e-2:
                        print cartographicMapArray[i][0]
                        listSurfacePoints = cartographicMapArray[i][2]
                        isInside = WindingNumber(p, listSurfacePoints)
                        if abs(isInside) > 5e-1:
                            #print 'is inside ', isInside
                            ini = rasterInsideCartographicArray[y][x]
                            if ini == 0:
                                rasterInsideCartographicArray[y][x] = [[cartographicMapArray[i][0],cartographicMapArray[i][1]]]
                            else: 
                                ini.append([cartographicMapArray[i][0],cartographicMapArray[i][1]])
                                rasterInsideCartographicArray[y][x] = ini
            #print rasterInsideCartographicArray[y][x]            
            break
        if (y == width_01) : 
            print '1%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
            break
        if (y == width_05) : 
            print '5%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_10) : 
            print '10%% donein %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_20) : 
            print '20%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_30) : 
            print '30%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_40) : 
            print '40%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_50) : 
            print '50%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_60) :
            print '60%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_70) : 
            print '70%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_80) : 
            print '80%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_90) : 
            print '90%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        if (y == width_100): 
            print '100%% done in %s seconds' % str(round(time.time()-t0,0))
            t0 = time.time()
        break    
    print '1 line done in %s seconds' % str(round(time.time()-t0,0))
    print('Saving the info cartographic into file')
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    fl = open(path_geo + "infoCartographic.csv", 'w')
    writer = csv.writer(fl, delimiter=';')
    for values in rasterInsideCartographicArray:
        writer.writerow(values)
    fl.close()       
    

""" Main Code """
#rasterToCoordinatesIni()
readInitialMap()
readRasterData(file = 'barcelona_raster_augusto_3000x3000 v9.asc')
readInitialMapRegionClassified(n_elements = 10)
#calculateRasterInsideCartographic2()



"""
# Create a jpeg file with the new interpolated velocity 
from osgeo import gdalnumeric
meshMeanVelArrayInteger = meshMeanVelArray.astype(gdalnumeric.uint8)
gdalnumeric.SaveArray(meshMeanVelArrayInteger, "/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/vel2.jpeg", format='JPEG')
m_x = meshPixelArray[0:10,0:10,0]
m_y = meshPixelArray[0:10,0:10,1]
m_x.tofile("/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/x_points.csv", sep=";")
m_y.tofile("/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/y_points.csv", sep=";")
meshMeanVelArray.tofile("/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/vel_points.csv", sep=";")
"""
