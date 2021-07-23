# -*- coding: utf-8 -*-
from keras.applications.vgg19 import VGG19
from keras.preprocessing import image
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten
import numpy as np
import concurrent.futures
import time
import tensorflow as tf
tf.compat.v1.enable_eager_execution()
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
import random
import openslide
from keras.utils import multi_gpu_model
import skimage.io
import pandas as pd
from collections import Counter
import warnings
import libtiff
from tqdm import tqdm
libtiff.libtiff_ctypes.suppress_warnings()
warnings.filterwarnings('ignore')
#====01_85_liver_images---------------
#clinical=pd.read_excel('/data/backup/Yuni/CRC_Liver/06_clinical_info/clinical_info_all.xlsx',sheet_name='first_85_liver')
#liver_name=list(clinical.iloc[:,0])
#f=open('/data/backup/Yuni/CRC_Liver/06_clinical_info/clinical_first_85_liver.pkl','wb')
#pickle.dump([liver_name],f)
#f.close()

f=open('/data/backup/Yuni/CRC_Liver/06_clinical_info/clinical_first_85_liver.pkl','rb')
liver_name=pickle.load(f)
f.close()
#liver_name=list(itertools.chain.from_iterable(liver_name))
liver_name=[str(name) for name in liver_name[0]]

os.chdir('/data/backup/Yuni/CRC_Liver/')
#os.environ['CUDA_VISIBLE_DEVICES'] =sys.argv[1]
os.environ['CUDA_VISIBLE_DEVICES'] ='3,4,5'
config = tf.compat.v1.ConfigProto()

#from tensorflow.python.client import device_lib
#
#def get_available_gpus():
#    local_device_protos = device_lib.list_local_devices()
#    return [x.name for x in local_device_protos if x.device_type == 'GPU']
#tf.test.is_gpu_available()

config.gpu_options.per_process_gpu_memory_fraction = 0.9
config.gpu_options.allow_growth = True
start = time.perf_counter()


#===================multi-gpu-model=======
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
parallel_model = multi_gpu_model(model, gpus=3)
parallel_model.load_weights('/data/backup/Yuni/CRC_Liver/04_fintune_train/model_finetuned/my_fitune_model_oct_05.h5')
parallel_model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])
parallel_model.summary()
#===========================================

labels = ['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK','Mucus']
#
#input_folder = '/data/backup/Yuni/CRC_Liver/liver_SVS/'
#tile_folder='/data/backup/Yuni/CRC_Liver/05_tiles_output/'
#output_folder='/data/backup/Yuni/CRC_Liver/07_new_output/'

input_folder = '/data/backup/Yuni/CRC_Liver/00_SVS_images/01_SVS_image_85/'
tile_folder='/data/backup/Yuni/CRC_Liver/01_tiles_output/01_tiles_output/'
output_folder='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_single_85/01_classified_first/'

#input_folder = '/data/backup/Yuni/CRC_Liver/08_new_SVS/'
#tile_folder='/data/backup/Yuni/CRC_Liver/08_1_tiles_output/'
#output_folder='/data/backup/Yuni/CRC_Liver/09_output_multi_liver/'

slidelist = os.listdir(tile_folder)
slidelist2=[]
for slide in slidelist:
    if os.path.exists(output_folder+slide.split('.')[0]+".pkl")==False:
        slidelist2.append(slide)
    else:
        print(f'{slide} exit!')

slidelist=slidelist2

slidelist2=[]
for name in slidelist:
    #print(name.split('.')[0])
    #name=slidelist[1]
    if name.split('.')[0] in liver_name:
        slidelist2.append(name)
    else:
        print(name)
slidelist=slidelist2

print("unclassify svs images number: ",len(slidelist))   
#print(int(sys.argv[2]),int(sys.argv[3]))         
#slidelist = slidelist[int(sys.argv[2]):int(sys.argv[3])]
#allslide_res=[]

#slide='210-2 (1).svs'

#"/data/backup/Yuni/CRC_Liver/07_new_output/210-2 (1).pkl"

def load_tiles(tile):
    tileP=f'{tile_folder}{slide}/{tile}'
    img = image.load_img(tileP, target_size=(224, 224))
    x = image.img_to_array(img)
    x = np.expand_dims(x, axis=0)
    #images.append(x)
    #print(f'Done loading...{tile}')
    return (x)


def run(load_tiles, tiles):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(load_tiles, tiles), total=len(tiles)))
    return results


def prediction(slideN):
    #slideN=slidelist[0]
    print(f'predicting {slideN}......')
    slide=slideN
    sliName=slide.split('.')[0]
    outPath =f'{output_folder}{sliName}'
    #if not os.path.exists(outPath):
        #os.makedirs(outPath)
#    slideImage = openslide.open_slide(f'{input_folder}{slide}')
    slide_res=[]
    tiles=os.listdir(f'{tile_folder}{slide}')
    #print(tiles)
    images=[]
    print(f'----loading {slideN}......')
    results=run(load_tiles, tiles)
    
#    with concurrent.futures.ThreadPoolExecutor() as executor:
#        results =executor.map(load_tiles, tiles)
    for res in results:
        images.append(res)
    images = np.vstack(images)
    print('---deep learning begin-----')
    classes = parallel_model.predict(images,batch_size=256,verbose=1)
    print('---deep learning finished-----')
    #print(classes)
    predicted_classes = np.argmax(classes, axis=1)
    name = [tile.split('.')[0] for tile in tiles]
    label=[]
    [label.append(labels[predicted_class]) for predicted_class in predicted_classes]
    slide_res.append([name,predicted_classes,label])
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
            #try:
            res=prediction(slide)
            print(res)
#            except:
#                print(f'{slide} failed!!!!!!!!!!!!!!!')
#                continue
        else:
            print(f'{slide} exit!')
