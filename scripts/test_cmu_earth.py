# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 15:01:37 2024

@author: sanghao
"""
import numpy as np
import pyproj
import scipy.spatial.transform     

def geodetic2enu(lat, lon, alt, lat_org, lon_org, alt_org):
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    x, y, z = transformer.transform( lon,lat,  alt,radians=False)
    print(x, y, z)
    x_org, y_org, z_org = transformer.transform( lon_org,lat_org,  alt_org,radians=False)
    vec=np.array([[ x-x_org, y-y_org, z-z_org]]).T
    
    print("delta ecef", vec)
    
    rot1 =  scipy.spatial.transform.Rotation.from_euler('x', -(90-lat_org), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rot3 =  scipy.spatial.transform.Rotation.from_euler('z', -(90+lon_org), degrees=True).as_matrix()#angle*-1 : left handed *-1

    rotMatrix = rot1.dot(rot3)    
    
    print(rotMatrix)
    
    enu = rotMatrix.dot(vec).T.ravel()
    return enu.T

def enu2geodetic(x,y,z, lat_org, lon_org, alt_org):
    transformer1 = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    transformer2 = pyproj.Transformer.from_crs(
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        )
    
    x_org, y_org, z_org = transformer1.transform( lon_org,lat_org,  alt_org,radians=False)
    ecef_org=np.array([[x_org,y_org,z_org]]).T
    
    rot1 =  scipy.spatial.transform.Rotation.from_euler('x', -(90-lat_org), degrees=True).as_matrix()#angle*-1 : left handed *-1
    rot3 =  scipy.spatial.transform.Rotation.from_euler('z', -(90+lon_org), degrees=True).as_matrix()#angle*-1 : left handed *-1

    rotMatrix = rot1.dot(rot3)

    ecefDelta = rotMatrix.T.dot( np.array([[x,y,z]]).T )
    ecef = ecefDelta+ecef_org
    
    print(ecef)
    
    lon, lat, alt = transformer2.transform( ecef[0,0],ecef[1,0],ecef[2,0],radians=False)

    return [lat,lon,alt]



if __name__ == '__main__':
    # The local coordinate origin (Zermatt, Switzerland)
    lat_org = 30 # deg
    lon_org = 120  # deg
    alt_org = 0     # meters

    # The point of interest
    lat = 30.001  # deg
    lon = 120.001   # deg
    alt = 10      # meters

    res1 = geodetic2enu(lat, lon, alt, lat_org, lon_org, alt_org)
    print (res1)

    
    x=10
    y=20
    z=30
    res2 = enu2geodetic(x,y,z, lat_org, lon_org, alt_org)
    print (res2)
