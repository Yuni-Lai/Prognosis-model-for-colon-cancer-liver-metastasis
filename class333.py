# -*- coding: utf-8 -*-
from keras.applications.vgg19 import VGG19
from keras.preprocessing import image
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten
import numpy as np
import concurrent.futures
import time
import tensorflow as tf
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
os.chdir('/data/backup/Yuni/CRC_Liver/')
os.environ['CUDA_VISIBLE_DEVICES'] = sys.argv[1]#'7'
config = tf.compat.v1.ConfigProto()

#from tensorflow.python.client import device_lib
#
#def get_available_gpus():
#    local_device_protos = device_lib.list_local_devices()
#    return [x.name for x in local_device_protos if x.device_type == 'GPU']
tf.test.is_gpu_available()

config.gpu_options.per_process_gpu_memory_fraction = 0.9
config.gpu_options.allow_growth = True
start = time.perf_counter()


base_model = VGG19(include_top=False, weights=None, input_shape=(224, 224, 3))
tmpx = base_model.output
tmpx = Flatten(name='flatten')(tmpx)
tmpx = Dense(4096,  activation='relu', name='fc1')(tmpx) # input_shape=(25088,),
tmpx = Dropout(0.5)(tmpx)
tmpx = Dense(4096,  activation='relu', name='fc2')(tmpx) # input_shape=(4096,),
tmpx = Dropout(0.5)(tmpx)
predictions = Dense(7,  activation='softmax')(tmpx)

model = Model(inputs=base_model.input, outputs=predictions)
model.load_weights('/data/backup/Yuni/CRC_Liver/04_fintune_train/model_finetuned/my_fitune_model_june26.h5')
# ??,??????????(???????) ???? vgg19 ????
for layer in base_model.layers:
    layer.trainable = True

# ????(??????????)
sgd = SGD(lr=1e-3, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])


#parallel_model = multi_gpu_model(model, gpus=1)
#parallel_model.compile(loss='categorical_crossentropy',optimizer=sgd, metrics=['accuracy'])
#parallel_model.summary()
labels = ['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK','Mucus']

input_folder = '/data/backup/Yuni/CRC_Liver/08_new_SVS/'
tile_folder='/data/backup/Yuni/CRC_Liver/08_1_tiles_output/'
output_folder='/data/backup/Yuni/CRC_Liver/09_output_multi_liver/'
slidelist = os.listdir(tile_folder)
#allslide_res=[]
slidelist2=slidelist[int(sys.argv[2]):int(sys.argv[3])]
#slide='210-2 (1).svs'

#"/data/backup/Yuni/CRC_Liver/07_new_output/210-2 (1).pkl"
def prediction(slideN):
    print(f'predicting {slideN}......')
    slide=slideN
    sliName=slide.split('.')[0]
    outPath =f'{output_folder}{sliName}'
    #if not os.path.exists(outPath):
        #os.makedirs(outPath)
    
    slideImage = openslide.open_slide(f'{input_folder}{slide}')
    slide_res=[]
    tiles=os.listdir(f'{tile_folder}{slide}')
    #print(tiles)
    images=[]
    for i,tile in enumerate(tiles):
        tileP=f'{tile_folder}{slide}/{tile}'
        img = image.load_img(tileP, target_size=(224, 224))
        x = image.img_to_array(img)
        x = np.expand_dims(x, axis=0)
        images.append(x)
    images = np.vstack(images)
    classes = model.predict(images,batch_size=30)
    #print(classes)
    predicted_classes = np.argmax(classes, axis=1)
    name = [tile.split('.')[0] for tile in tiles]
    label=[]
    [label.append(labels[predicted_class]) for predicted_class in predicted_classes]
    slide_res.append([name,predicted_classes,label])
    f=open(f'{outPath}.pkl','wb')
    pickle.dump([name,predicted_classes,label],f)
    f.close()
#    f=open(f'{outPath}.pkl','rb')
#    name,predicted_classes,label=pickle.load(f)
#    f.close()
    #slideImage = openslide.open_slide('D:/00_Apr05_CRC_Liver_Image/SVS_image/424498.svs')

    dim=slideImage.level_dimensions[1]
    a1=int(dim[0]/224)
    a2=int(dim[1]/224) 

    #simg = slideImage.get_thumbnail((a1,a2))

    #plt.figure()
    #plt.imshow(simg)
    maxxx=max(a1,a2)
    index_max=np.argmax(dim)
    classmap=np.zeros([maxxx,maxxx])
    x=np.zeros(len(name))
    y=np.zeros(len(name))
    for i,n in enumerate(name):#'0_0'
        a,b=n.split('_')
        a,b=int(a),int(b)
        classmap[b,a]=predicted_classes[i]+1#['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK']=1~9; 0=white back ground
        x[i],y[i]=b,a
    if index_max==0:      
        classmap=classmap[0:a2,:]
    else:
        classmap=classmap[:,0:a1]
#    f=open(f'{outPath}_classmap_x_y.pkl','wb')
#    pickle.dump([classmap,x,y],f)
#    f.close()
#    f=open(f'{outPath}_classmap_x_y.pkl','rb')
#    classmap,x,y=pickle.load(f)
#    f.close() 
    #label = ['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK','Mucus']
   
    #colors =['red','purple','brown','green','blue','gray']
    #colors =[(105,105,105), (255,255,255), (114,64,70), (195,100,197), (252,108,133), (205,92,92), (255,163,67), (70,130,180), (28,172,120)]
    #simg=np.array(simg)
    #plt.imshow()
    #classmap_rgb = color.label2rgb(classmap,image=simg,colors=colors,image_alpha=0, alpha=1)#,bg_label=0,bg_color=bg_color
    #classmap_rgb = color.label2rgb(classmap,colors=colors)#,bg_label=0,bg_color=bg_color
    #classmap_rgb = color.label2rgb(classmap)
    #plt.style.use('classic')
    #plt.figure()
    #plt.imshow(classmap_rgb)
    #skimage.io.imsave(f'{outPath}_res.jpg', classmap_rgb)
#    x=np.array([0,1,2,3,4,5,6,7,8,9])
#    x=np.repeat(x, 3)
#    x=np.tile(x,(3,1))
#    legend = color.label2rgb(x,colors=colors,alpha=1,image_alpha=0,kind='overlay')
#    plt.figure()
#    plt.imshow(legend)
#    plt.xlim(-1,30)
#    plt.ylim(-10,30)
    #plt.text(0.3, -1, 'ADI', fontdict={'size': 8, 'color': 'r'})
#    label2 = ['BACK','ADI','BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']
#    for i,text in enumerate(label2):
#        plt.text(0.3+i*3, -1.5, f'{text}', fontdict={'size': 12, 'color': 'r'}) 
#    STR_Tum_slide=[]
#    LYM_Tum_slide=[]
#    NORM_Tum_slide=[]
#    MUS_Tum_slide=[]
#    MUC_Tum_slide=[]
#    ADI_Tum_slide=[]
#    DEB_Tum_slide=[]
#    
#    a,b=classmap.shape
#    #------only sample from Tumor areas----------
##    Tumor_mask=(classmap==9)
#    Tumor_index=np.where(classmap==1)
#    #plt.imshow(Tumor_mask)
#    random.seed(10)
#    sample_N=len(Tumor_index[1])
#    num=list(range(0,sample_N))
#    random_num = random.sample(num, sample_N)
#    Sample_num=int(sample_N*0.5)
#    Tumor_index=np.array(Tumor_index)  
#    random_index=Tumor_index[:,random_num[1:Sample_num]]
#    #STR_mask=(classmap==8)
#    #plt.imshow(STR_mask)
#    for index in random_index.transpose().tolist():
#        i=index[0]
#        j=index[1]
#        #print(i,j)
#        if i-5>0 and i+5<a and j-5>0 and j+5<b:
#            submap=classmap[i-5:i+5,j-5:j+5]
#            #print(i,j)
#            #print(submap)
#            submap=submap.tolist()
#            submap = list(itertools.chain(*submap))
#            counts=np.zeros(7)
#            for k in range(1,7):#1~9,no'0'
#                #print(k)
#                counts[k] = submap.count(k)#['BACK==0','ADI','BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']
#                #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK']
#            del submap
#            #print(counts)
#            All=sum(counts)
#            All_interest=All-counts[2]
#            if counts[1]/All>=0.3:#
#                STR_Tum=counts[2]/All_interest
#                STR_Tum_slide.append(STR_Tum)
#                #print(STR_Tum_slide)
#                LYM_Tum=counts[4]/All_interest
#                LYM_Tum_slide.append(LYM_Tum)
#                NORM_Tum=counts[5]/All_interest
#                NORM_Tum_slide.append(NORM_Tum)
#                DEB_Tum=counts[3]/All_interest
#                DEB_Tum_slide.append(DEB_Tum)
#    All_proportion_slide=[STR_Tum_slide,LYM_Tum_slide,DEB_Tum_slide,NORM_Tum_slide]
#    Tum_score=[np.mean(pro) for pro in All_proportion_slide]
    #============radar graph==========================

#    plt.rcParams['font.sans-serif'] = 'Microsoft YaHei'
#    plt.rcParams['axes.unicode_minus'] = False
#    
#
#    plt.style.use('ggplot')
#    
#
#    values = Tum_score
#    feature = ['STR_TUM','LYM_TUM','MUS_TUM','MUC_TUM','ADI_TUM','DEB_TUM','NORM_TUM']
#    
#    N = len(values)
#
#    angles=np.linspace(0, 2*np.pi, N, endpoint=False)
#    
#
#    values=np.concatenate((values,[values[0]]))
#    angles=np.concatenate((angles,[angles[0]]))
#    
#
#    fig=plt.figure()
#
#    ax = fig.add_subplot(111, polar=True)
#
#    ax.plot(angles, values, 'o-', linewidth=2)
#
#    ax.fill(angles, values, alpha=0.25)
#
#    ax.set_thetagrids(angles * 180/np.pi, feature)
#
#    ax.set_ylim(0,1)
#
#    plt.title('Tumor interection')
#
#    ax.grid(True)
#    plt.savefig(f'{outPath}_interaction.jpg')

    #plt.show()
    #====================================================
    f=open(f'{outPath}.pkl','wb')
    pickle.dump([name,predicted_classes,label,classmap,x,y,classes],f)
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
    for slide in slidelist2:
        if os.path.exists(output_folder+slide.split('.')[0]+".pkl")==False:
            try:
                res=prediction(slide)
                print(res)
            except:
                print(f'{slide} failed!!!!!!!!!!!!!!!')
                continue
        
    