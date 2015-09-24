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
from PIL import Image
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import numpy as np
import os.path
import math
from osgeo.gdalconst import *


# Find the different information (translation, rotation and scalefactor) needed to change pixels to coordinates 
def pixelsToCoordinatesIni(P1,P2):
    global transMatrix, transVector, scaleFactor, originPointPixel, originPointCoord
    #Initial coordinates of 2 points (lower and lefter points: 1 and 2)
    print 'Calculating the needed information (translation, rotation and scalefactor) to change coordinates-raster to pixel-velocity'
    p1 = {  'x': 129515.0,
            'y': -29683.0   }
    p2 = {  'x': 135770.0,
            'y': -14263.0   } 
            
#    p1 = {  'x': 928271.0,
#            'y': 4587398.0  }
#    p2 = {  'x': 933208.0,
#            'y': 4603698.2  } 
            


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
    

# Reading the cartographic map. The info is stored into cartographicMapArray(parcelId,RiskueId,ListCoorXY,center)
def readInitialMap(filename='Predios_bcn_ETRS89.srl'):
    global cartographicMapArray
    cartographicMapArray = []
    path_cad = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Cadastro/'
    if os.path.isfile(path_cad + filename):
        print 'Reading serialized cadastre from file %s' % filename
        file = open(path_cad + filename, 'r')
        i = 0
        for line in file:
            i += 1
            #Seeking the info related to coordinates that define the polygon
            if line.find(';') >= 0:
                coordinates = line.split(';')[-1]
                coordinates = coordinates.split('|')[0:-1]
                #print 'coor: ', coordinates
                newCoor = [[float(coordinates[n]), float(coordinates[n+1])] for n in xrange(0,len(coordinates),2)]
                #coor = [[float(coordinates[n]), float(coordinates[n+1])] for n in xrange(0,len(coordinates),2)]
                #newCoor = readInitialMapCheckPoints(coor)
                #if len(coor) > len(newCoor): print 'removed from parcel: %s %s of %s' % (parcelId,len(newCoor),len(coor))
                #print 'coor new: ', newCoor
                center = np.mean(newCoor, axis=0)
                #print 'center: ', center
            else:
                parcelId = line.split('|')[4]
                parcelId = parcelId.replace('#RISKUE', '')
                RiskueId = line.rstrip().split('|')[-1]
                #print 'parcel: ', parcelId
                #print 'riskue: ', RiskueId
                cartographicMapArray.append([parcelId,RiskueId,newCoor,center.tolist()])
            #if i > 3: break
    else: 
        print 'file %s does not exist' % str(path_cad + filename)      


#Checking the list of point and return only the needed points (to avoid points closer to another)
def readInitialMapCheckPoints(listPoints):
    newListPoints = listPoints
    for i in listPoints:
        for j in listPoints:
            print i, ' ', j, ' ', calculateDistance(i,j)
            if i != j and calculateDistance(i,j) < 2.0:
                print 'removed point: %s' % j
                newListPoints.remove(j)
    return newListPoints
    

# Reading geometry data from the raster file (values are stored in a matrix called geoArray)
def readRasterData(file = 'barcelona_raster_augusto_500x400.asc'):
    global geoArray, top_left_x, resolution_x, top_left_y, resolution_y
    #initial variables
    path_geo = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Geo Raster/'
    ini_height = 386   # Value needed because the raster file is displaced in z
    if os.path.isfile(path_geo + file):
        print 'Reading initial geometry from file %s' % file
        src_ds = gdal.Open(path_geo + file)
        srcband = src_ds.GetRasterBand(1)  
        cols = src_ds.RasterXSize
        rows = src_ds.RasterYSize
        NoDataValue = srcband.GetNoDataValue()
        geoTransform = src_ds.GetGeoTransform()
        top_left_x = geoTransform[0] # top left x 
        resolution_x = geoTransform[1] # w-e resolution 
        top_left_y = geoTransform[3] # top left y
        resolution_y = geoTransform[5] # n-s resolution 
        geoArray = srcband.ReadAsArray(0, 0, cols, rows)
        # Considering the initial_height
        geoArray = geoArray - ini_height
        # Changing NoDataValue into dataArray to nan
        geoArray[geoArray == NoDataValue] = nan
    else:
        print 'file %s does not exist' % str(path_geo + file)



# Calculating the solid angle of a point respect to one surface. If the solid angle is 0, that means that the
# point is outside, if the value is 2PI the point is inside, and any other value the point is over one of the 
# edges of the surface.
def WindingNumber (p, ListSurfacePoints):
    totalAngle = 0
    """The points are normalized using p as origin amb scaled using 1/1000 """
    #qx = p[0]
    #qy = p[1] 
    qx = 0.0
    qy = 0.0 
    for i in range(1,len(ListSurfacePoints)):
        p1 = ListSurfacePoints[i-1]
        p2 = ListSurfacePoints[i]
        #cx_p1 = p1[0]
        #cy_p1 = p1[1]
        #cx_p2 = p2[0]
        #cy_p2 = p2[1] 
        cx_p1 = (p1[0]-p[0])*0.001
        cy_p1 = (p1[1]-p[1])*0.001
        cx_p2 = (p2[0]-p[0])*0.001
        cy_p2 = (p2[1]-p[1])*0.001  
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


#Return the distance between two points
def calculateDistance(p, centralPoint):
    if not isinstance(p, ndarray):
        p = np.array(p)
    if not isinstance(centralPoint, ndarray):
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
            """tranformating coordinates is needed"""
            row_list.append((x, y))
        rasterCoordArray.append(row_list) 

    rasterCoordArray = np.array(rasterCoordArray)    

    print('Calculating the raster points inside cartographic')
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
    for y in np.arange(dim_Y):
        for x in np.arange(dim_X):
            print x, ' ', y
            # Find the edges of the square in Coord and
            p = rasterCoordArray[x][y]
            print p
            #then transform it into UTM coordinates
            #p = coordinatesToPixels(p)   
            p = [932912.591045, 4601664.363319]
            #and then I look for the inside points and points over the edges using the windingNumber
            for i in range(0,len(cartographicMapArray)):
                center = cartographicMapArray[i][3]
                distance = calculateDistance(center, p)
                if distance <= minDistance:
                    print cartographicMapArray[i][0]
                    listSurfacePoints = cartographicMapArray[i][2]
                    isInside = WindingNumber(p, listSurfacePoints)
                    if abs(isInside) > 5e-1:
                        print 'is inside ', isInside
                        ini = rasterInsideCartographicArray[x][y]
                        if ini == 0:
                            rasterInsideCartographicArray[x][y] = [[cartographicMapArray[i][0],cartographicMapArray[i][1]]]
                        else: 
                            ini.append([cartographicMapArray[i][0],cartographicMapArray[i][1]])
                            rasterInsideCartographicArray[x][y] = ini
            print ini            
            break
        break
        if (y == width_10) : print '10% done'
        if (y == width_20) : print '20% done'
        if (y == width_30) : print '30% done'
        if (y == width_40) : print '40% done'
        if (y == width_50) : print '50% done'
        if (y == width_60) : print '60% done'
        if (y == width_70) : print '70% done'
        if (y == width_80) : print '80% done'
        if (y == width_90) : print '90% done'
        if (y == width_100): print '100% done'



def writeVelRaster(file='0-1.tiff'):
    path = '/home/daniel/Documentos/Ofertes/Recurs Eolic/Estudi/Dades inicials/Mapes minieolica/'
    path_file = path + file
    
    # My array lon / lat
    dim_X, dim_Y, dim_Z = meshCoordArray.shape
    latArray = np.empty((dim_X,dim_Y), dtype='f')
    lonArray = np.empty((dim_X,dim_Y), dtype='f')    
    for y in range(dim_Y):
        for x in range(dim_X):
            val = meshCoordArray[x][y].tolist()
            latArray[x][y] = val[1]
            lonArray[x][y] = val[0]
   
    # For each pixel I know it's latitude and longitude.
    # As you'll see below you only really need the coordinates of
    # one corner, and the resolution of the file.
    
    xmin,ymin,xmax,ymax = [lonArray.min(),latArray.min(),lonArray.max(),latArray.max()]
    nrows,ncols = dim_X, dim_Y
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform=(xmin,xres,0,ymax,0, -yres)   
    # That's (top left x, w-e pixel resolution, rotation (0 if North is up), 
    #         top left y, rotation (0 if North is up), n-s pixel resolution)
    # I don't know why rotation is in twice???
    
    output_raster = gdal.GetDriverByName('GTiff').Create(path_file,ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)                     # This one specifies WGS84 lat long.
                                                 # Anyone know how to specify the 
                                                 # IAU2000:49900 Mars encoding?
    output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system 
                                                       # to the file
    output_raster.GetRasterBand(1).WriteArray(meshMeanVelArray)
        
    

""" Main Code """
readInitialMap()
readRasterData()
calculateRasterInsideCartographic()


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



