# -*- coding: utf-8 -*-
from keras.applications.vgg19 import VGG19
from keras.preprocessing import image
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten
import numpy as np
#import tensorflow as tf

import os
import pickle
from skimage import color
import matplotlib.pyplot as plt
import itertools
#from statistics import median

import random
import openslide

import skimage.io
os.chdir('/data/backup/Yuni/CRC_Liver/')
#os.environ['CUDA_VISIBLE_DEVICES'] = '0,1'
#config = tf.ConfigProto()
#config.gpu_options.per_process_gpu_memory_fraction = 0.9

base_model = VGG19(include_top=False, weights=None, input_shape=(224, 224, 3))

tmpx = base_model.output
tmpx = Flatten(name='flatten')(tmpx)
tmpx = Dense(4096,  activation='relu', name='fc1')(tmpx)
tmpx = Dropout(0.5)(tmpx)
tmpx = Dense(4096,  activation='relu', name='fc2')(tmpx)
tmpx = Dropout(0.5)(tmpx)
predictions = Dense(9,  activation='softmax')(tmpx)

model = Model(inputs=base_model.input, outputs=predictions)
model.load_weights('./9Structure_bestweight.hd5')  # 9 structures

labels = ['ADI', 'BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']

input_folder = '/data/backup/Yuni/CRC_Liver/SVS_image/'
tile_folder='/data/backup/Yuni/CRC_Liver/01_tiles_output/'
output_folder='/data/backup/Yuni/CRC_Liver/02_classified_output/'
slidelist = os.listdir(tile_folder)
#allslide_res=[]
#slide='424498.svs'


def prediction(slideN):
    slide=slideN
    sliName=slide.split('.')[0]
    outPath =f'{output_folder}{sliName}'
    #if not os.path.exists(outPath):
        #os.makedirs(outPath)
    print('===========1========')
    slideImage = openslide.open_slide(f'{input_folder}{slide}')
    slide_res=[]
    tiles=os.listdir(f'{tile_folder}{slide}')
    #print(tiles)
    images=[]
    print('===========2========')
    for i,tile in enumerate(tiles):
        tileP=f'{tile_folder}{slide}/{tile}'
        print(tileP)
        img = image.load_img(tileP, target_size=(224, 224))
        x = image.img_to_array(img)
        x = np.expand_dims(x, axis=0)
        images.append(x)
    print('===========3========')
    images = np.vstack(images)
    classes = model.predict(images,batch_size=100)
    #print(classes)
    predicted_classes = np.argmax(classes, axis=1)
    name = [tile.split('.')[0] for tile in tiles]
    label=[]
    [label.append(labels[predicted_class]) for predicted_class in predicted_classes]
    slide_res.append([name,predicted_classes,label])
    print('===========4========')
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
    print('===========5========')
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
        classmap[b,a]=predicted_classes[i]+1#['ADI', 'BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']=1~9; 0=white back ground
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
    label = ['ADI', 'BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']
   
    colors =['linen','purple','gray','brown','green','olive','cyan','blue','yellow', 'red']
    #colors =[(105,105,105), (255,255,255), (114,64,70), (195,100,197), (252,108,133), (205,92,92), (255,163,67), (70,130,180), (28,172,120)]
    #simg=np.array(simg)
    #plt.imshow()
    #classmap_rgb = color.label2rgb(classmap,image=simg,colors=colors,image_alpha=0, alpha=1)#,bg_label=0,bg_color=bg_color
    classmap_rgb = color.label2rgb(classmap,colors=colors)#,bg_label=0,bg_color=bg_color
    #classmap_rgb = color.label2rgb(classmap)
    #plt.style.use('classic')
    #plt.figure()
    #plt.imshow(classmap_rgb)
    print('===========6========')
    skimage.io.imsave(f'{outPath}_res.jpg', classmap_rgb)
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
    STR_Tum_slide=[]
    LYM_Tum_slide=[]
    NORM_Tum_slide=[]
    MUS_Tum_slide=[]
    MUC_Tum_slide=[]
    ADI_Tum_slide=[]
    DEB_Tum_slide=[]
    
    a,b=classmap.shape
    #------only sample from Tumor areas----------
#    Tumor_mask=(classmap==9)
    Tumor_index=np.where(classmap==9)
    #plt.imshow(Tumor_mask)
    random.seed(10)
    sample_N=len(Tumor_index[1])
    num=list(range(0,sample_N))
    random_num = random.sample(num, sample_N)
    Sample_num=int(sample_N*0.5)
    Tumor_index=np.array(Tumor_index)  
    random_index=Tumor_index[:,random_num[1:Sample_num]]
    #STR_mask=(classmap==8)
    #plt.imshow(STR_mask)
    for index in random_index.transpose().tolist():
        i=index[0]
        j=index[1]
        #print(i,j)
        if i-5>0 and i+5<a and j-5>0 and j+5<b:
            submap=classmap[i-5:i+5,j-5:j+5]
            #print(i,j)
            #print(submap)
            submap=submap.tolist()
            submap = list(itertools.chain(*submap))
            counts=np.zeros(10)
            for k in range(1,10):#1~9,no'0'
                #print(k)
                counts[k] = submap.count(k)#['BACK==0','ADI','BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']
            del submap
            #print(counts)
            All=sum(counts)
            All_interest=All-counts[2]
            if counts[9]/All>=0.5:#
                STR_Tum=counts[8]/All_interest
                STR_Tum_slide.append(STR_Tum)
                #print(STR_Tum_slide)
                LYM_Tum=counts[4]/All_interest
                LYM_Tum_slide.append(LYM_Tum)
                NORM_Tum=counts[7]/All_interest
                NORM_Tum_slide.append(NORM_Tum)
                MUS_Tum=counts[6]/All_interest
                MUS_Tum_slide.append(MUS_Tum)
                MUC_Tum=counts[5]/All_interest
                MUC_Tum_slide.append(MUC_Tum)
                ADI_Tum=counts[1]/All_interest
                ADI_Tum_slide.append(ADI_Tum)
                DEB_Tum=counts[3]/All_interest
                DEB_Tum_slide.append(DEB_Tum)
    All_proportion_slide=[STR_Tum_slide,LYM_Tum_slide,MUS_Tum_slide,MUC_Tum_slide,ADI_Tum_slide,DEB_Tum_slide,NORM_Tum_slide]
    Tum_score=[np.mean(pro) for pro in All_proportion_slide]
    #============radar graph==========================

    plt.rcParams['font.sans-serif'] = 'Microsoft YaHei'
    plt.rcParams['axes.unicode_minus'] = False
    

    plt.style.use('ggplot')
    

    values = Tum_score
    feature = ['STR_TUM','LYM_TUM','MUS_TUM','MUC_TUM','ADI_TUM','DEB_TUM','NORM_TUM']
    
    N = len(values)

    angles=np.linspace(0, 2*np.pi, N, endpoint=False)
    

    values=np.concatenate((values,[values[0]]))
    angles=np.concatenate((angles,[angles[0]]))
    

    fig=plt.figure()

    ax = fig.add_subplot(111, polar=True)

    ax.plot(angles, values, 'o-', linewidth=2)

    ax.fill(angles, values, alpha=0.25)

    ax.set_thetagrids(angles * 180/np.pi, feature)

    ax.set_ylim(0,1)

    plt.title('Tumor interection')

    ax.grid(True)
    print('===========7========')
    plt.savefig(f'{outPath}_interaction.jpg')

    #plt.show()
    #====================================================
    f=open(f'{outPath}.pkl','wb')
    pickle.dump([name,predicted_classes,label,classmap,x,y,classmap_rgb,Tum_score],f)
    f.close()
    return f'Done Classification...{slideN}'


#slide='424498.svs'
#prediction(slide)
from concurrent.futures import ProcessPoolExecutor
import time
if __name__ == "__main__":
  start = time.perf_counter()
  p = ProcessPoolExecutor(max_workers=40)
  results = p.map(prediction, slidelist)
  p.shutdown(wait=True)
  for result in results:
      print(result)
  
  finish = time.perf_counter()
  
  print(f'Finished in {round(finish-start, 2)} second(s)')   
    