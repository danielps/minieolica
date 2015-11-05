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
import os.path
import math
import csv
import time

 

# Reading the cartographic map. The info is stored into cartographicMapArray(parcelId,RiskueId,ListCoorXY,center)
def readInitialMap(filename='Predios_bcn_ETRS89.srl'):
    global cartographicMapArray
    t0 = time.time()
    cartographicMapArray = []
    parcelsToNotBeConsidered = [5398001,6173015,6210005,6760019,6770005,6770007,6770009,6780009,
                                6780015,8670006,8690006,11000007,11040007,11040008,16850001,
                                17050001,20022001,20321001,20342001,20342002,20343001,20803003,20804002,
                                20804004,20804006,20833002,21581001,22113001,22710001,22729023,22771001,
                                22780002,22919013,22919013,28421018,28583017,28611001,28611002,29642001,
                                29694004,29694009,29694010,29722004,29722004,29722055,29723016,33656001,
                                33736075,34568077,34569006,34586003,34586004,34625003,34673017,34681007,
                                34682036,34686011,34686064,34725020,34815001,34824058,34824061,35600017,
                                35611022,35611026,40154011,60860001,81604001,83962018,84042030,84406020,
                                90771001,90773002,94815058,96500003,96500039,96500059,96500073,96500085,
                                96500086,96500095,97050005,97071001,97071002,97071006,97071026,97081001,
                                97082004,97294001,97394001,97406002,98100006,98110003,98121025,98180001,99535001]

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
                if RiskueId != '0' and parcelId not in parcelsToNotBeConsidered:
                    newCoor = readInitialMapCheckPoints(newCoor)
                    cartographicMapArray.append([parcelId,RiskueId,newCoor,center.tolist(),square,dist/2.0])
            #if i > 3: break
    else: 
        print 'file %s does not exist' % str(path_cad + filename)      
    
    print 'Reading serialized cadastre has finished in %s seconds' % str(round(time.time()-t0,1))


# Reading the cartographic map. The info is stored into cartographicMapArray(parcelId,RiskueId,ListCoorXY,center)
def readInitialReducedMap(filename='Predios_bcn_ETRS89.srl', listParcels =[]):
    t0 = time.time()
    cartographicReducedMapArray = []
 
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
                if parcelId in listParcels:
                    newCoor = readInitialMapCheckPoints(newCoor)
                    cartographicReducedMapArray.append([parcelId,RiskueId,newCoor,center.tolist(),square,dist/2.0])
            #if i > 3: break
    else: 
        print 'file %s does not exist' % str(path_cad + filename)      
    return cartographicReducedMapArray
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
    
#Returns the list of index (cartographicMapRegionIndex) to a given point
def findMapRegion(point):
    dim_X , dim_Y = mapRegion.shape 
    for y in np.arange(dim_Y):
        for x in np.arange(dim_X):
            square = mapRegion[x][y]
            isInsideSquare = WindingNumber(point, square)
            if abs(isInsideSquare) > 5e-2:
                return cartographicMapRegionIndex[x][y]   
    


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
    intervalPercentage = 5
    listPercentage = [int(dim_Y*x/100.0) for x in range(0,101,intervalPercentage)]
    t0 = time.time()
    for y in np.arange(dim_Y):
        for x in np.arange(dim_X):
            #print x, ' ', y
            # Find the edges of the square in Coord and
            p = rasterCoordArray[y][x]
            #print p
            #Look for the mapRegion where the point belongs to 
            mapRegionListIndex = findMapRegion(p)
            #and then I look for the inside points and points over the edges using the windingNumber
            for i in mapRegionListIndex:
                center = cartographicMapArray[i][3]
                distance = calculateDistance(center, p)
                minDistance = cartographicMapArray[i][5]
                if distance <= minDistance:
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
        if y in listPercentage:
            print '%s%% done in %s seconds' % (str(y*100/dim_Y), str(round(time.time()-t0,0)))
            t0 = time.time()
        #break    
    print('Saving the info cartographic into file')
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    fl = open(path_geo + "infoCartographic.csv", 'w')
    fl.write('NCOLS %s \n' %dim_X)
    fl.write('NROWS %s \n' %dim_Y)
    fl.write('XLLCORNER %s \n' %top_left_x)
    fl.write('YLLCORNER %s \n' %(top_left_y+resolution_y*dim_Y))
    fl.write('CELLSIZE %s \n' %resolution_x)
    fl.write('NODATA_VALUE 0 \n')
    writer = csv.writer(fl, delimiter=';')
    for values in rasterInsideCartographicArray:
        writer.writerow(values)
    fl.close()       
    
    
#Update the info with the new parcels to be considered
def updateRasterInsideCartographic(infoParcels):
    print('Calculating the coordinates of the raster file')  
    dim_X , dim_Y = geoArray.shape
    updateRasterInsideCartographicArray = np.zeros((dim_X,dim_Y), dtype=object) 
    rasterCoordArray = []
    for y in np.arange(top_left_y, top_left_y+resolution_y*dim_Y, resolution_y):
        row_list = []
        for x in np.arange(top_left_x, top_left_x+resolution_x*dim_X, resolution_x):
            row_list.append((x, y))
        rasterCoordArray.append(row_list) 

    rasterCoordArray = np.array(rasterCoordArray)    
    print('Calculating the raster points inside cartographic')
    intervalPercentage = 5
    listPercentage = [int(dim_Y*x/100.0) for x in range(0,101,intervalPercentage)]
    t0 = time.time()
    for y in np.arange(dim_Y):
        for x in np.arange(dim_X):
            #print x, ' ', y
            # Find the edges of the square in Coord and
            p = rasterCoordArray[y][x]
            #print p
            #and then I look for the inside points and points over the edges using the windingNumber
            for i in infoParcels:
                center = i[3]
                distance = calculateDistance(center, p)
                minDistance = i[5]
                if distance <= minDistance:
                    #print cartographicMapArray[i][0]
                    listSurfacePoints = i[2]
                    isInside = WindingNumber(p, listSurfacePoints)
                    if abs(isInside) > 5e-1:
                        #print 'is inside ', isInside
                        ini = updateRasterInsideCartographicArray[y][x]
                        if ini == 0:
                            #print i[0]
                            updateRasterInsideCartographicArray[y][x] = [[i[0],i[1]]]
                        else: 
                            ini.append([i[0],i[1]])
                            updateRasterInsideCartographicArray[y][x] = ini
            #print rasterInsideCartographicArray[y][x]  
        if y in listPercentage:
            print '%s%% done in %s seconds' % (str(y*100/dim_Y), str(round(time.time()-t0,0)))
            t0 = time.time()
        #break    
    print('Saving the info cartographic into file')
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    fl = open(path_geo + "infoCartographic_new.csv", 'w')
    writer = csv.writer(fl, delimiter=';')
    for values in updateRasterInsideCartographicArray:
        writer.writerow(values)
    fl.close()       
    
    

def updateCSV():
    from numpy import genfromtxt
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    my_data = genfromtxt(path_geo + "infoCartographic_new.csv", delimiter=';',dtype=object)
    dim_X , dim_Y = geoArray.shape
    my_new_data = np.zeros((dim_X,dim_Y), dtype=object) 
    for y in range(0,3000):
        for x in range(0,3000):
            if my_data[y][x] != '0':
                print y, x, my_data[y][x]
    my_ori_data = genfromtxt(path_geo + "infoCartographic.csv", delimiter=';',dtype=object)
    for y in range(0,3000):
        for x in range(0,3000):
            if my_data[y][x] == '0' and my_ori_data[y][x] != '0':
                my_new_data[y][x] = my_ori_data[y][x]
            elif my_data[y][x] != '0' and my_ori_data[y][x] == '0':
                my_new_data[y][x] = my_data[y][x]
            elif my_data[y][x] != '0' and my_ori_data[y][x] != '0':
                my_new_data[y][x] = my_data[y][x]                
                #print 'error ', my_data[y][x], my_ori_data[y][x]
            else:
                my_new_data[y][x] = '0'
    print('Saving the new info cartographic into file')
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    fl = open(path_geo + "infoCartographic_new_def.csv", 'w')
    writer = csv.writer(fl, delimiter=';')
    for values in my_new_data:
        writer.writerow(values)
    fl.close()       
 

""" Main Code """

#readInitialMap()
#readRasterData(file = 'barcelona_raster_augusto_3000x3000 v9.asc')
#readInitialMapRegionClassified(n_elements = 10)
#calculateRasterInsideCartographic()




def prova(filename = 'raster_prova_new.asc'):

    parcelsToNotBeConsidered = ['05070010','05398001','05995001','05996001','06066005','06069001','06173015','06210005',
                                '06240016','06251004','06280012','06280017','06280019','06280020','06721024','06760003',
                                '06760019','06770001','06770002','06770003','06770004','06770005','06770007',
                                '06770008','06770009','06780008','06780009','06780010','06780011','06780015','06780016',
                                '06780017','07013019','07091015','07110004','07127003','07162005','07162007','07222011',
                                '07234010','07234019','07234015','08083016','08112039','08131023','08220002','08220005','08412007',
                                '08413003','08413004','08413005','08413011','08580001','08670006','08680003','08680004',
                                '08680005','08680006','08690003','08690006','08750001','08760001','08920001','08940001',
                                '08960001','08970001','10780011','10780014','10790001','11000007','11040007','11040008',
                                '11900001','11920004','11920006','16541025','16850001','17050001','20090004','20110002',
                                '20022001','20321001','20342001','20342002','20343001','20803003','20804002','20804004',
                                '20804006','20833002','21570001','21581001','21602001','21602004','21615002','21830003',
                                '21830016','21830018','21840004','21870012','21902003','21902015','22091015','22113001',
                                '22302015','22210029','22701001','22710001','22729001','22729023',
                                '22771001','22780002','22904008','22918001','22919001','22919013','23090078','23110017',
                                '27293009','27293030','27293081','27293277','27335001','27971012','28005013','28010008',
                                '28250001','28370013','28421018','28421019','28553002','28581018','28583017','28583018',
                                '28583019','28583020','28590003','28590012','28590013','28590014','28590017','28611001',
                                '28611002','28612001','28612003','28612005','28612012','28626003','28632008','28632010',
                                '28632012','28643005','28643012','28656004','28656005','28656008','28657001','28659008',
                                '29642001','29694004','29694009','29694010','29697001','29722004','29722055','29723016',
                                '32230001','33656001','33736075','33990022','34130018','34141013','34568077','34569006',
                                '34586003','34586004','34625003','34630001','34630004','34630006','34673017','34681007',
                                '34682036','34686011','34686064','34725020','34815001','34824058','34824061','35600017',
                                '35611022','35611026','40154011','60860001','70433025','71723034','82000007','81604001',
                                '83594035','83594057','83962018','84042030','84406020','90771001','90773002','96330005',
                                '94815058','96500003','96500039','96500059','96500073','96500085','96500086','96500095',
                                '97050005','97071001','97071002','97071006','97071026','97081001','97082004','97294001',
                                '97394001','97406002','98100006','98110003','98121025','98180001','98180008','99535001']
                        

                                
    parcelsToBeConsidered = ['08040021','08050002','08050003','08131020','08131021','08450007','08450015','08450016',
                             '08450017','08450019','08450020','08450021','08450024','08450026','08450027','08420021',
                             '08450028','08420022','08460014','08433001','08433002','08433003','08433004','08433005',
                             '08433006','08433007','08433009','08433010','08433011','08433012','08433013','08433014',
                             '08433015','08475017','08433018','08433019','08433020','08433021','08433022','08433023',
                             '08433024','08433025','08433026','08433027','08433029','08433030','08433031','08433032',
                             '08433033','08433035','08433036','08433037','08433038','08433039','08433040','08433041',
                             '08433042','08433043','08433044','08433045','08433046','08433047','08433048','08433049',
                             '08433050','08443001','08443002','08443003','08443004','08443005','08443006','08443009',
                             '08443010','08443012','08443013','08443015','08443016','08443017','08443018','08443019',
                             '08443020','08443021','08443022','08443023','08443024','08443025','08443026','08443027',
                             '08443028','08443029','08443030','08443031','08443032','08443033','08443034','08443035',
                             '08443036','08443037','08443038','08443039','08443040','08443041','08443042','08443043',
                             '08443044','08443045','08443046','08443047','08443048','08443050','08561009','08580004',
                             '22324002','22336056','22336046','22336052','22336049','22336047','22336020','22336018',
                             '22336019','22336016','22336017','22336015','22336014','22336011','22336012','28626006',
                             '22918014','23020015','72744048','72744049','72744050','72744053','72744002','72744054',
                             '72744056','72744071','72744072','72744073','72744074','72744035','72744037','72744038',
                             '72744039','72744043','72744044','72744045','72744046','72744047',]                             
    parcelsToBeConsidered = []
    #   '06912001','06911001','07070007','07091009','07091010','07091006','07091008','07091005','07776003'
                             
    if len(parcelsToBeConsidered) > 0:
        readRasterData(file = 'barcelona_raster_augusto_3000x3000 v9.asc')
        infoParcels = readInitialReducedMap(filename='Predios_bcn_ETRS89.srl', listParcels = parcelsToBeConsidered)
        updateRasterInsideCartographic(infoParcels)
        updateCSV()
        
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    #path_geo = '/home/xcipriano/Escriptori/MiniEolica/Geo Raster/'
    fileread =  "infoCartographic_new_def.csv"
    fr = open(path_geo + fileread, 'r')
    fl = open(path_geo + filename, 'w')
    fl.write('NCOLS 3000 \n')
    fl.write('NROWS 3000 \n')
    fl.write('XLLCORNER %s \n' %top_left_x)
    fl.write('YLLCORNER %s \n' %(top_left_y+resolution_y*3000))
    fl.write('CELLSIZE %s \n' %resolution_x)
    fl.write('NODATA_VALUE 9999 \n')
    writer = csv.writer(fl, delimiter=' ')
    for line in fr:
        transf_line = []
        new_line = line.split(';')
        #print len(new_line)
        for elem in new_line:
            if elem == '0' or elem == '0\r\n':
                transf_line.append('9999')
            else: 
                ok = []
                #print elem
                elem_split = elem.split("'")
                new_elem_list = []
                new_elem_list.append(elem_split[1])
                if len(elem_split) > 5:
                    new_elem_list.append(elem_split[5])
                if len(elem_split) > 9: 
                    new_elem_list.append(elem_split[9])
                if len(elem_split) > 13: 
                    new_elem_list.append(elem_split[13])        
                if len(elem_split) > 17:
                    print elem_split
                for i in new_elem_list:
                    #print i
                    #I only eliminate the parcelelsnottobeconsidered if there is no other parcel that contains the point
                    if i in parcelsToNotBeConsidered:
                        ok.append(False)
                        #print 'False', i
                    else:
                        ok.append(True)                        
                if True in ok:  
                    #print 'i', elem_split
                    transf_line.append(elem_split[1])
                else:
                    transf_line.append('9999')
                #print len(transf_line)
        writer.writerow(transf_line)
        
    fr.close()
    fl.close() 
    
    
    
def fixFinalFile(filename = 'raster_prova_new.asc'):
    #Reading the original file
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    if os.path.isfile(path_geo + filename):
        print 'Reading initial parcel from file %s' % filename
        from numpy import genfromtxt
        fr = open(path_geo + filename, 'r')
        l = {}
        for i in range(6):
            l[i] = fr.readline()
        NoDataValue = int(l[5].split(' ')[1])
        parcelArray = genfromtxt(path_geo + filename, delimiter=' ',skip_header=6)
        parcelArray[parcelArray == NoDataValue] = np.nan
        dim_X, dim_Y = parcelArray.shape
        for x in range(1,dim_X-1):
            for y in range(1,dim_Y-1):
                if parcelArray[x,y] is np.nan:
                    reducedParcelArray = parcelArray[x-1:x+1,y-1:y+1]
    else:
        print 'Reading initial parcel, file %s does not exist' % str(path_geo + filename)