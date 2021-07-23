# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 12:42:24 2020

@author: Yuni
"""
from keras.applications.vgg19 import VGG19
from keras.preprocessing import image
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten
import numpy as np
import concurrent.futures
import time
import tensorflow as tf
tf.compat.v1.disable_eager_execution()
tf.compat.v1.disable_v2_behavior()
import sys
from keras.optimizers import SGD
import os
import pickle
from skimage import color
import matplotlib.pyplot as plt
import itertools
#from statistics import median
import concurrent.futures
from os import listdir
import random
import openslide
from keras.utils import multi_gpu_model
import skimage.io
import pandas as pd
from collections import Counter
from os.path import isfile, join
from openslide.deepzoom import DeepZoomGenerator
from normalizeStaining_my import normalizeStaining_my
import warnings
import libtiff
libtiff.libtiff_ctypes.suppress_warnings()
warnings.filterwarnings('ignore')
#clinical=pd.read_excel('/data/backup/Yuni/CRC_Liver/06_clinical_info/clinical_info_all.xlsx',sheet_name='liver')

#liver_name=list(clinical.iloc[:,0])
os.chdir('/data/backup/Yuni/CRC_Liver/')
os.environ['CUDA_VISIBLE_DEVICES'] = '0,1,2,3,4'
config = tf.compat.v1.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.9

tf.test.is_gpu_available()

config.gpu_options.per_process_gpu_memory_fraction = 0.9
config.gpu_options.allow_growth = True
start = time.perf_counter()

#===================multi-gpu-model
base_model = VGG19(include_top=False, weights=None, input_shape=(224, 224, 3))
tmpx = base_model.output
tmpx = Flatten(name='flatten')(tmpx)
tmpx = Dense(4096,  activation='relu', name='fc1')(tmpx) # input_shape=(25088,),
tmpx = Dropout(0.5)(tmpx)
tmpx = Dense(4096,  activation='relu', name='fc2')(tmpx) # input_shape=(4096,),
tmpx = Dropout(0.5)(tmpx)
predictions = Dense(7,  activation='softmax')(tmpx)
sgd = SGD(lr=1e-3, decay=1e-6, momentum=0.9, nesterov=True)
with tf.device("/cpu:0"):
    model = Model(inputs=base_model.input, outputs=predictions)

for layer in base_model.layers:
    layer.trainable = True

parallel_model = multi_gpu_model(model, gpus=5)
parallel_model.load_weights('/data/backup/Yuni/CRC_Liver/04_fintune_train/model_finetuned/my_fitune_model_oct_05.h5')
parallel_model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])
parallel_model.summary()
#===================



labels = ['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK','Mucus']
#
#input_folder = '/data/backup/Yuni/CRC_Liver/liver_SVS/'
#tile_folder='/data/backup/Yuni/CRC_Liver/05_tiles_output/'
#output_folder='/data/backup/Yuni/CRC_Liver/07_new_output/'

input_folder = '/data/backup/Yuni/CRC_Liver/SVS_image/'
#tile_folder='/data/backup/Yuni/CRC_Liver/08_1_tiles_output/'
output_folder='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_classified_first/'

image_names=os.listdir(input_folder)
image_names2=[]
for slide in image_names:
    
    if os.path.exists(output_folder+slide.split('.')[0]+".pkl")==False:
        image_names2.append(slide)
    else:
        print(f'{slide} exit!')

slidelist=image_names2


print("unclassify svs images number: ",len(slidelist))   
#print(int(sys.argv[2]),int(sys.argv[3]))         
#slidelist = slidelist2[int(sys.argv[2]):int(sys.argv[3])]
#allslide_res=[]

#slide='210-2 (1).svs'

#"/data/backup/Yuni/CRC_Liver/07_new_output/210-2 (1).pkl"
pixel = 224 # or 256
pixel_num = int(pixel*pixel*0.5) # *3/2

def warp(tile_loc):
    i,j=tile_loc
    tile = data_gen.get_tile(num_level - 2, (i, j)).convert("RGB")
    img_np_1 = np.array(tile)
    if img_np_1.shape == (pixel, pixel, 3):
        IBW = np.sum([img_np_1[:, :, 0], img_np_1[:, :, 1], img_np_1[:, :, 2]], axis=0) / 3
        if np.sum(IBW >= 225) < pixel_num and np.sum(IBW <= 20) < pixel_num * 1.5:
            #print(i,j)
            tile_temp=normalizeStaining_my(img_np_1, savePath=None)
            a,b,_=tile_temp.shape
            if a==224 and b==224:
                tile_temp = image.img_to_array(tile_temp)
                tile_temp = np.expand_dims(tile_temp, axis=0)
                images.append(tile_temp)
                tilesN.append([i,j])
                return f'Done normalizing...{i,j}'
            
def prediction(slideN):
    #slideN=slidelist[0]
    print(f'predicting {slideN}......')
    slide=slideN
    sliName=slide.split('.')[0]
    outPath =f'{output_folder}{sliName}'
    #if not os.path.exists(outPath):
        #os.makedirs(outPath)
    
#    slideImage = openslide.open_slide(f'{input_folder}{slide}')
    image_path = os.path.join(input_folder, slideN)
    slide = openslide.open_slide(image_path)
    data_gen = DeepZoomGenerator(slide, tile_size=pixel, overlap=0, limit_bounds=True)
    num_level = data_gen.level_count
    num_tiles = data_gen.level_tiles[num_level - 2]  # range from 0  # [num_level - 2] equals to 20x
    #data_gen.level_dimensions[num_level - 2]        # [num_level-1] equals to 40x
    x_steps = num_tiles[0]  # number of x starting points
    y_steps = num_tiles[1]  # number of y starting points
    images=[]
    tilesN=[]
    tiles_locations=[]
    for i in range(x_steps):
        for j in range(y_steps):
            tiles_locations.append([i,j])
            
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(warp, tiles_locations)
        
    for i in range(x_steps):
        for j in range(y_steps):
            # get a patch
            #i,j=10,180
            tile = data_gen.get_tile(num_level - 2, (i, j)).convert("RGB")
            img_np_1 = np.array(tile)
            if img_np_1.shape == (pixel, pixel, 3):
                IBW = np.sum([img_np_1[:, :, 0], img_np_1[:, :, 1], img_np_1[:, :, 2]], axis=0) / 3
                if np.sum(IBW >= 225) < pixel_num and np.sum(IBW <= 20) < pixel_num * 1.5:
                    #print(i,j)
                    try:
                        tile_temp=normalizeStaining_my(img_np_1, savePath=None)
                        a,b,_=tile_temp.shape
                        if a==224 and b==224:
                            x = image.img_to_array(tile_temp)
                            x = np.expand_dims(x, axis=0)
                            images.append(x)
                            tilesN.append([i,j])
                    except:
                        continue    
    
    images = np.vstack(images)
    classes = model.predict(images,batch_size=64)
    #print(classes)
    predicted_classes = np.argmax(classes, axis=1)
    name = tilesN
    label=[]
    [label.append(labels[predicted_class]) for predicted_class in predicted_classes]
    #slide_res.append([name,predicted_classes,label])
    f=open(f'{outPath}.pkl','wb')
    pickle.dump([name,predicted_classes,label],f)
    f.close()
    return f'Done Classification...{slideN}'


#if __name__ == "__main__":
#    #results = executor.map(do_something, secs)
#    p = ProcessPoolExecutor(max_workers=5)
#    results = p.map(do_something, secs)
#    p.shutdown(wait=True)
#    for result in results:
#        print(result)
#    
#    finish = time.perf_counter()
#    
#    print(f'Finished in {round(finish-start, 2)} second(s)')
    
#if __name__ == "__main__":
#    print(slidelist)
#    with concurrent.futures.ProcessPoolExecutor() as executor:
#        results = executor.map(prediction, slidelist)
##    p = ProcessPoolExecutor(max_workers=3)
##    results = p.map(prediction, slidelist)
##    p.shutdown(wait=True)
#    for result in results:
#        print(result)
#    
#    finish = time.perf_counter()
#    
#    print(f'Finished in {round(finish-start, 2)} second(s)')   
#    
    
if __name__ == "__main__":
    for slide in slidelist:
        if os.path.exists(output_folder+slide.split('.')[0]+".pkl")==False:
            try:
                res=prediction(slide)
                print(res)
            except:
                print(f'{slide} failed!!!!!!!!!!!!!!!')
                continue
        else:
            print(f'{slide} exit!')
        
    