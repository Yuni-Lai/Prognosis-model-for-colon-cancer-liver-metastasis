# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 09:37:03 2020
@author: Yuni
"""
import os
import cv2
import pickle
import numpy as np
from skimage import color
import matplotlib.pyplot as plt
import pandas as pd
from skimage import transform
from skimage import measure
from skimage import morphology
import copy
import json
import openslide
import libtiff
import ipdb
libtiff.libtiff_ctypes.suppress_warnings()
#with open('D:/00_Apr05_CRC_Liver_Image/python/template_CRC2.json') as f:
#    data = json.load(f)
#print(data)
#f=open('D:/00_Apr05_CRC_Liver_Image/python/data_template.txt','wb')
#pickle.dump(data,f)
#f.close()
f=open('/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/data_template.txt','rb')
data=pickle.load(f)
f.close()
tile_output_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/05_beijing/01_classified/'
#'/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_classified_first/'
svs_path='/data/backup2/Yuni/CRC_Liver/00_SVS_images/05_beij_511/2020-10-15/'
#'/data/backup/Yuni/CRC_Liver/SVS_image/'
qupath_json_output='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/05_beijing/02_qupath_res/'
#'/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/02_qupath_res/'
reslist = os.listdir(tile_output_path)
def get_json_annotation_by_classmap(res):
    
    #res='424498.pkl'
    #res='425876.pkl'
    f=open(f'{tile_output_path}{res}','rb')
    name,predicted_classes,label=pickle.load(f)
    f.close()
    slideImage = openslide.open_slide(f"{svs_path}{res.split('.')[0]}.svs")
    dim=slideImage.level_dimensions[1]
    a1=int(dim[0]/224)*2+1
    a2=int(dim[1]/224)*2+1
    maxxx=max(a1,a2)
    index_max=np.argmax(dim)
    classmap=np.zeros([maxxx,maxxx])
    x=np.zeros(len(name))
    y=np.zeros(len(name))
    for i,n in enumerate(name):#'0_0'
        a,b=n[0],n[1]
        a,b=int(a),int(b)
        classmap[b,a]=predicted_classes[i]+1#['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK']=1~9; 0=white back ground
        x[i],y[i]=b,a
    if index_max==0:      
        classmap=classmap[0:a2,:]
    else:
        classmap=classmap[:,0:a1]
    #Tumscore_all.append(Tum_score)

    classmap2=np.array(classmap)
    #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK'==6,'Mucus'==7]
    #colors =['linen','red','yellow','brown','green','cyan','linen','blue']
    #colors =['linen','purple','gray','brown','green','olive','cyan','blue','yellow', 'red']
    #colors =[(0.41,0.41,0.41), (1,1,1), (0.33,0.25,0.27), (0.76470588, 0.39215686, 0.77254902),(0.98823529, 0.42352941, 0.52156863), (0.80392157, 0.36078431, 0.36078431), (1, 0.63921569, 0.2627451), (0.2745098 , 0.50980392, 0.70588235), (0.10980392, 0.6745098 , 0.47058824)]#,
    #classmap_rgb = color.label2rgb(classmap2,colors=colors,kind='overlay')
    #    plt.figure()
    #plt.imshow(classmap_rgb)
    #ipdb.set_trace()
    typeids=[1,2,3,4,5,7]
    json_list=[]
    for typeid in typeids:
        print(typeid)
        #typeid=1
        #ipdb.set_trace()
        Mask=(classmap2==typeid)
    #        plt.imshow(Mask)
        #Mask=morphology.remove_small_holes(Mask, 5)
    #        plt.figure()
    #        plt.imshow(Mask1)
        Mask=morphology.remove_small_objects(Mask,5)
    #        plt.figure()
    #        plt.imshow(Mask2)
        #contours = measure.find_contours(Mask2, 0.5)
        Mask=Mask.astype('uint8')
        #cv2.RETR_EXTERNAL#cv2.RETR_TREE
        contours, hierarchy=cv2.findContours(Mask,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)#cv2.RETR_TREE
        #img=np.zeros((a2,a1))
        #cv2.drawContours(img, contours, -1, (255,255,255), 1)
        #plt.imshow(img)
#        contour=None
#        for con in contours:
#            #print(con)
#            if contour is None:
#                contour = con
#            else:
#                contour=np.vstack((contour,con))
        
        #contour=[contour]
#    #    show the contour with mask
#        fig, ax = plt.subplots()
#        ax.imshow(Mask2,cmap=plt.cm.gray)
#        
#        for n, contour in enumerate(contours):
#            ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
#        
#        ax.axis('image')
#        ax.set_xticks([])
#        ax.set_yticks([])
#        plt.show()
    #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK']
        try:
            if typeid==1:
                for contour in contours:
                    #i=0
                    #contour=contours[i]
                    contour=np.squeeze(contour,axis = 1)
                    if contour.shape[0]>=4:
                        contour= contour*448
                        contour= contour.astype(np.int)
                        #contour[:,[0,1]] = contour[:,[1,0]]
                        #contour[:,]=contour[:,]+
                        contour= contour
                        temp=copy.deepcopy(data[3])
                        #ipdb.set_trace()
                        temp['geometry']['coordinates'].append(contour.tolist())
                        json_list.append(temp)
            if typeid==2:
                for contour in contours:
                    contour=np.squeeze(contour,axis = 1)
                    if contour.shape[0]>=4:
                        contour= contour*448
                        contour= contour.astype(np.int)
                        #contour[:,[0,1]] = contour[:,[1,0]]
                        contour= contour
                        temp=copy.deepcopy(data[2])
                        temp['geometry']['coordinates'].append(contour.tolist())
                        json_list.append(temp)
            if typeid==3:
                for contour in contours:
                    contour=np.squeeze(contour,axis = 1)
                    if contour.shape[0]>=4:
                        contour= contour*448
                        contour= contour.astype(np.int)
                        #contour[:,[0,1]] = contour[:,[1,0]]
                        contour= contour
                        temp=copy.deepcopy(data[1])
                        temp['geometry']['coordinates'].append(contour.tolist())
                        json_list.append(temp)
            if typeid==4:
                for contour in contours:
                    contour=np.squeeze(contour,axis = 1)
                    if contour.shape[0]>=4:
                        contour= contour*448
                        contour= contour.astype(np.int)
                        #contour[:,[0,1]] = contour[:,[1,0]]
                        contour= contour
                        temp=copy.deepcopy(data[0])
                        temp['geometry']['coordinates'].append(contour.tolist())
                        json_list.append(temp)
            if typeid==5:
                for contour in contours:
                    contour=np.squeeze(contour,axis = 1)
                    if contour.shape[0]>=4:
                        contour= contour*448
                        contour= contour.astype(np.int)
                        #contour[:,[0,1]] = contour[:,[1,0]]
                        contour= contour
                        temp=copy.deepcopy(data[5])
                        temp['geometry']['coordinates'].append(contour.tolist())
                        json_list.append(temp)
            if typeid==7:
                for contour in contours:
                    contour=np.squeeze(contour,axis = 1)
                    if contour.shape[0]>=4:
                        contour= contour*448
                        contour= contour.astype(np.int)
                        #contour[:,[0,1]] = contour[:,[1,0]]
                        contour= contour
                        temp=copy.deepcopy(data[4])
                        temp['geometry']['coordinates'].append(contour.tolist())
                        json_list.append(temp)
        except:
            continue
    #        for con in contour:
    #            [str(int(con[0]*448)), str(int(con[1]*448))]
    for json1 in json_list:
        if json1['geometry']['coordinates'][0][0] != json1['geometry']['coordinates'][0][-1]:
            json1['geometry']['coordinates'][0].append(json1['geometry']['coordinates'][0][0])
    num=res.split(".")[0]
    with open(f'{qupath_json_output}{num}.json', 'w') as f:
        json.dump(json_list,f)
        print(f"Generate .Json annotation of {res} finished!")
#    a,b=classmap.shape
#    dst=transform.resize(classmap, (a*224, b*224))
#    plt.imshow(dst)
    
for res in reslist[0:10]:
    print(res)
    get_json_annotation_by_classmap(res)