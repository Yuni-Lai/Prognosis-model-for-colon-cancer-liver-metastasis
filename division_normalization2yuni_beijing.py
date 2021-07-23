import openslide
import os
import numpy as np
import sys
from keras.preprocessing import image
from normalizeStaining_my import normalizeStaining_my
from threading import Thread
from openslide.deepzoom import DeepZoomGenerator
from os import listdir
from os.path import isfile, join
import pandas as pd
from collections import Counter
import concurrent.futures
from tqdm import tqdm
import libtiff
import pickle
libtiff.libtiff_ctypes.suppress_warnings()
os.environ["OMP_NUM_THREADS"] = '8'
#clinical=pd.read_excel('/data/backup/Yuni/CRC_Liver/06_clinical_info/clinical_info_all.xlsx',sheet_name='liver')

#liver_name=list(clinical.iloc[:,0])

#sys.path.append('D:/00_Apr05_CRC_Liver_Image')
input_folder = '/data/backup2/Yuni/CRC_Liver/00_SVS_images/05_beij_511/2020-10-15/'
output_folder = '/data/backup/Yuni/CRC_Liver/01_tiles_output/05_tiles_output/'
output_folder2 = '/data/backup2/Yuni/CRC_Liver/01_tiles_output/05_tiles_output/'
output_folder3='/dataserver145/image/Yuni/CRC_Liver/01_tiles_output/05_tiles_output/'
#
#just get the name of images in a folder
#def get_image_names(folder):
#	image_names = [f for f in listdir(folder) if isfile(join(folder, f))]
#	if '.DS_Store' in image_names:
#		image_names.remove('.DS_Store')
#	image_names = sorted(image_names)
#	return image_names

def get_image_names(folder):
    imagefolder=os.listdir(input_folder)
    image_names=[]
    for f in imagefolder:
        images=os.listdir(input_folder+f)
        images = [f+'/'+im for im in images]
        image_names.extend(images)
    return image_names


image_names = os.listdir(input_folder)
#tile_names = os.listdir(input_folder)

#image_names_liver=[]
#for name in image_names:
#    #print(name.split('.')[0])
#    if name.split('.')[0] in liver_name:
#        image_names_liver.append(name)
#    else:
#        print(name)
#image_names=image_names_liver

#C1=Counter(liver_name)#265
#C2=Counter(image_names_liver)#262



pixel = 224 # or 256
pixel_num = int(pixel*pixel*0.5) # *3/2


def cut_tiles(tile_loc):
    #tile_loc=tiles_locations[0]
    i,j=tile_loc
    tile = data_gen.get_tile(num_level - 2, (i, j)).convert("RGB")
    img_np_1 = np.array(tile)
    if img_np_1.shape == (pixel, pixel, 3):
        IBW = np.sum([img_np_1[:, :, 0], img_np_1[:, :, 1], img_np_1[:, :, 2]], axis=0) / 3
        if np.sum(IBW >= 225) < pixel_num and np.sum(IBW <= 20) < pixel_num * 1.5:
            try:
                tile_temp=normalizeStaining_my(img_np_1, savePath=None)
                a,b,_=tile_temp.shape
                if a==224 and b==224:
                    tile_temp = image.img_to_array(tile_temp)
                    tile_temp = np.expand_dims(tile_temp, axis=0)
                    #images.append(tile_temp)
                    #tilesN.append([i,j])
                    return (tile_loc,tile_temp)
            except:
                print("exception")
                pass
    
def run(cut_tiles, tiles_locations):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(cut_tiles, tiles_locations), total=len(tiles_locations),mininterval=30,maxinterval=60))
        #results = executor.map(cut_tiles, tiles_locations)
    return results


def wrapper(image_name):
    #image_name=image_names[0]
    #image_name='2020-06-29/744897_1A'
    Name=image_name.split('.')[0].replace(" ","_")
    #subfolder=Name.split('/')[0]
    if os.path.exists(f'{output_folder}{Name}.pkl')==True:
        return f'{Name} exit!'
    if os.path.exists(f'{output_folder2}{Name}.pkl')==True:
        return f'{Name} exit!'
    image_path = os.path.join(input_folder, image_name)
    slide = openslide.open_slide(image_path)
    #output_folder1 = (output_folder + image_name[k])
    output_folder1=output_folder3
    #if os.path.exists(output_folder1):
     #   print('file exit: ',output_folder1)
    
    if not os.path.exists(output_folder1):
        os.makedirs(output_folder1)
    global data_gen
    data_gen = DeepZoomGenerator(slide, tile_size=pixel, overlap=0, limit_bounds=True)
    
    global num_level
    num_level = data_gen.level_count
    
    num_tiles = data_gen.level_tiles[num_level - 2]  # range from 0  # [num_level - 2] equals to 20x
    #data_gen.level_dimensions[num_level - 2]        # [num_level-1] equals to 40x
    x_steps = num_tiles[0]  # number of x starting points
    y_steps = num_tiles[1]  # number of y starting points
    tiles_locations=[]
    for i in range(x_steps):
        for j in range(y_steps):
            tiles_locations.append([i,j])
    print(f'----loading {image_name}......')
    results=run(cut_tiles, tiles_locations)
    images=[]
    tilesN=[]
    for res in results:
        if res is not None:
            images.append(res[1])
            tilesN.append(res[0])
    f=open(f'{output_folder1}{Name}.pkl','wb')
    pickle.dump([images,tilesN],f)
    f.close()
#        
#        for i in range(x_steps):
#            for j in range(y_steps):
#                # get a patch
#                tile = data_gen.get_tile(num_level - 2, (i, j)).convert("RGB")
#                img_np_1 = np.array(tile)
#    
#                if img_np_1.shape == (pixel, pixel, 3):
#                    outPath = (output_folder1 + '/' + str(i) + "_" + str(j) + '.jpg')
#    
#                    IBW = np.sum([img_np_1[:, :, 0], img_np_1[:, :, 1], img_np_1[:, :, 2]], axis=0) / 3
#                    if np.sum(IBW >= 225) < pixel_num and np.sum(IBW <= 20) < pixel_num * 1.5:
#                        try:
#                            normalizeStaining_my(img_np_1, savePath=outPath)
#                        except:
#                            continue
    print("finished extraction and division: " + image_name)


#threads = []
#for i in range(5):
#    #argg=int(5*i)
#    t = Thread(target=wrapper, args=[i])
#    t.start()
#    threads.append(t)
#
#for thread in threads:
#    thread.join()

# parallel
#threads = []
#step=len(image_names)/5
#step_int=int(step)
#for i in range(step_int): 
#    t=Thread(target=wrapper, args=(5*i,5*i+5))
#    print(5*i,5*i+5)
#    threads.append(t)
#    t.start()
#
#t1=Thread(target=wrapper, args=(step_int*5,len(image_names)))
#print(step_int*5,len(image_names))
#threads.append(t1)
#t1.start()
#
#for thread in threads:
#    thread.join()


from concurrent.futures import ProcessPoolExecutor
import time

start = time.perf_counter()
if __name__ == "__main__":
    #results = executor.map(do_something, secs)
    p = ProcessPoolExecutor(max_workers=1)
    results = p.map(wrapper, image_names)
    p.shutdown(wait=True)
    for result in results:
        print(result)
    
    finish = time.perf_counter()
    
    print(f'Finished in {round(finish-start, 2)} second(s)')

