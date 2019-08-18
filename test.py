import numpy as np
import matplotlib.pyplot as plt
import rasterio
from OptimizedPitRemoval import PitRemoval

# global variable
RASTER= rasterio.open('test_dem.tif')
DATA= RASTER.read(1)

# test IsBorder function
def test_IsBorder():
    #read a rasterio map
    PR= PitRemoval('test_dem.tif')
    border_map= PR._init_queue(PR.dem).reshape(100,100)
    plt.figure()
    plt.imshow(border_map.astype(np.int32))

    return None

def test_Neighbors():
    PR= PitRemoval('test_dem.tif')
    neighbors= PR._get_neighbors(15)
    mat= np.zeros((3,3))
    mat[1,1]=15
    mat[0,:]=neighbors[:3]
    mat[1,2]=neighbors[3]
    mat[2,:]= neighbors[6:3:-1]
    mat[1,0]=neighbors[7]

    print(mat)

    return None

def test_main_queue():
    PR= PitRemoval('test_dem.tif')
    PR._init_queue(PR.dem)
    print(PR.geo_info['flow_dir'])
    print(len(PR.main_queue['id']))

def test_main_program():
    PR= PitRemoval('test_dem.tif')
    PR._init_queue(PR.dem)
    PR.iterate_main_queue()
    PR._write_dem(dst='cleaned.tif')

if __name__=='__main__':
    # test_IsBorder()   #pass
    # test_Neighbors().   #pass
    # test_main_queue().  #pass
    test_main_program()

