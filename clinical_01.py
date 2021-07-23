#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 18:29:14 2020

@author: Yuni
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 18:09:10 2020

@author: Yuni
"""


#import xlrd
import os
import pickle
import numpy as np
from skimage import color
import matplotlib.pyplot as plt
#import random
import pandas as pd
import itertools
import openslide
from skimage import morphology
from skimage import transform,data
import cv2
from scipy import stats
from skimage import measure
from sklearn import metrics
import seaborn as sns
from scipy import ndimage
from collections import Counter
import libtiff
libtiff.libtiff_ctypes.suppress_warnings()


#clinical=pd.read_excel('/data/backup/Yuni/CRC_Liver/06_clinical_info/clinical_info_all.xlsx',sheet_name='399_liver_survival')
##
#f=open('/data/backup/Yuni/CRC_Liver/06_clinical_info/399_liver_survival.pkl','wb')
#pickle.dump([clinical],f)
#f.close()


f=open('/data/backup/Yuni/CRC_Liver/06_clinical_info/399_liver_survival.pkl','rb')
clinical=pickle.load(f)
f.close()
clinical=clinical[0]
liver_name=[str(name) for name in clinical.iloc[:,1]]

#--------setting---path----------------
#01
global path_res
global path_svs
global root_path
root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_single_85/'
path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_single_85/01_classified_first/'
path_svs='/data/backup2/Yuni/CRC_Liver/00_SVS_images/01_SVS_image_85/'

#02
#global path_res
#global path_svs
#global root_path
#root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/02_single_285/'
#path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/02_single_285/01_classified_second/'
#path_svs='/data/backup2/Yuni/CRC_Liver/00_SVS_images/02_liver_SVS_285/'

#03
#global path_res
#global path_svs
#global root_path
#root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/03_mult_265/'
#path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/03_mult_265/01_classified_third/'
#path_svs='/data/backup2/Yuni/CRC_Liver/00_SVS_images/03_multi_regions_SVS_265/'

#04
#global path_res
#global path_svs
#global root_path
#root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/04_mult_141/'
#path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/04_mult_141/01_classified/'
#path_svs='/data/backup/Yuni/CRC_Liver/00_SVS_images/04_multi_regions/CRC_LM_10/'

#---rename svs files----
#slidelist = os.listdir(path_svs)
#newname=[slidename.replace(' ', '_') for slidename in slidelist]
#C1=Counter(newname)#578
#
#fail=[]
#for slidename in slidelist:
#    try:
#        os.rename(path_svs+slidename,path_svs+slidename.replace(' ', '_'))
#    except:
#        print(f'{slidename} failed!')
#        fail.append(slidename)
#        continue
#---------------------------------------------
#reslist1_new =[res.replace(' ', '') for res in reslist1]
#for i in range(len(reslist1)):
#    os.rename(path1+reslist1[i], path1+reslist1_new[i])

def get_image_names(input_folder):
    imagefolder=os.listdir(input_folder)
    image_names=[]
    for f in imagefolder:
        images=os.listdir(input_folder+f)
        images = [f+'/'+im for im in images]
        image_names.extend(images)
    return image_names

reslist = os.listdir(path_res)
#reslist = get_image_names(path_res)

reslist_new_number=[res.split('.pkl')[0] for res in reslist]

reslist_effect=[]
for res in reslist_new_number:
    if res in liver_name:
        reslist_effect.append(res)
    else:
        print(f'could not find {res}---')

#[name for name in liver_name if '63378' in name]
        
#Tumscore_all=[]
#classmap_rgb_all=[]

#import skimage

#res_patients=[res.split('_')[0] for res in reslist]
#patients=dict(Counter(res_patients))

#['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]
def get_tumor_size(res):
    try:
        f=open(f'{path_res}{res}.pkl','rb')
        name,predicted_classes,label=pickle.load(f)
        f.close()
        #slideImage = openslide.open_slide(f'{path_svs}{res}.svs')
  
        if os.path.getsize(f'{path_res}{res}.pkl') == 0:
            print(f'{path_res}{res}.pkl is empty, remove it!')
            os.remove(f'{path_res}{res}.pkl')
    except:
        print(f'something wrong with {path_res}{res}.pkl')
        return 0
    
    Tumor_mask=(predicted_classes==0)
    return(np.sum(Tumor_mask))



def calculate_feature2(res):#according to tumor region
#    try:
    #res='765680_7C'
    LYM_dis=[]
    STR_dis=[]
    print('processing:',res)
        #res=reslist1_new[0]
    try:
        f=open(f'{path_res}{res}.pkl','rb')
        #name,predicted_classes,label,classmap,x,y,classes=pickle.load(f)
        name,predicted_classes,label=pickle.load(f)
        f.close()
    except:
        print("file has something worng")
        return None
    slideImage = openslide.open_slide(f'{path_svs}{res}.svs')
    dim=slideImage.level_dimensions[1]#1
    a1=int(dim[0]/224)*2+1
    a2=int(dim[1]/224)*2+1
    maxxx=max(a1,a2)
    index_max=np.argmax(dim)
    classmap=np.zeros([maxxx,maxxx])
    x=np.zeros(len(name))
    y=np.zeros(len(name))
    for i,n in enumerate(name):#'0_0'
        a,b=n.split('_')
        #a,b=n[0],n[1]
        a,b=int(a),int(b)
        classmap[b,a]=predicted_classes[i]+1#['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK']=1~9; 0=white back ground
        x[i],y[i]=b,a
    if index_max==0:
        classmap=classmap[0:a2,:]
    else:
        classmap=classmap[:,0:a1]
    
    a,b=classmap.shape
#    colors=['linen','red','yellow','brown','green','cyan','linen','blue']
#    classmap_rgb = color.label2rgb(classmap,colors=colors,kind='overlay')
#    plt.imshow(classmap_rgb)
    #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]
    #寻找肿瘤区域，并计算STR和LYM到肿瘤区域的距离分布
    Tumor_mask=(classmap==1)
    #plt.imshow(Tumor_mask)
    #Tumor_mask=Tumor_mask.astype(np.uint8)*255
    kernel = morphology.disk(3)
    #kernel2 = np.ones((5,5),np.uint8) 
    Tumor_mask=morphology.remove_small_objects(Tumor_mask,50)
    Tumor_mask=Tumor_mask.astype(np.uint8)*255
    Tumor_mask = cv2.morphologyEx(Tumor_mask, cv2.MORPH_CLOSE, kernel)
    Tumor_mask = cv2.morphologyEx(Tumor_mask, cv2.MORPH_CLOSE, kernel)
    #Tumor_mask =morphology.closing(Tumor_mask, selem=None, out=None)
    #Tumor_mask=morphology.remove_small_holes(Tumor_mask, 1000)
    Tumor_mask=ndimage.binary_fill_holes(Tumor_mask)
    #Tumor_mask=(Tumor_mask==255)
    Tumor_mask=morphology.remove_small_objects(Tumor_mask,50)
    #plt.imshow(Tumor_mask)
    Tumor_mask=Tumor_mask.astype(np.uint8)*1

    contours, hierarchy = cv2.findContours(Tumor_mask, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)#Tumor area contours
    classes_for_analysis=['Stroma','Necrosis','Lymphocyte','Normal_hepatic','Mucus']#'Normal_hepatic'
    for class_name in classes_for_analysis:
        locals()[class_name+"_dis"]=[]
    #Applying polygon test to calculate the distance:
    class_number_relation_dict={'BACK':0,'Tumor':1,'Stroma':2,'Necrosis':3,'Lymphocyte':4,'Normal_hepatic':5,'BACK2':6,'Mucus':7}
    reversed_dict = {v : k for k, v in class_number_relation_dict.items()}
    reversed_dict[0]
    try:
        for i in range(a):
            for j in range(b):
                c=classmap[i,j]
                if c!=0 and c!=6 and c!=1:
                    temp=[]
                    for cnt in contours:
                        temp.append(cv2.pointPolygonTest(cnt, (i, j), True))
                    dis=max(temp)
                    locals()[reversed_dict[c]+"_dis"].append(dis)
    except:
        print(f'{res} failed!')
        return None
    LYM_dis.append(eval(reversed_dict[4]+"_dis"))
    STR_dis.append(eval(reversed_dict[2]+"_dis"))
    #print(reversed_dict[4])
    classmap_list=classmap.tolist()
    classmap_list = list(itertools.chain(*classmap_list))
    counts_all=np.zeros(8)
    for k in range(1,8):#1~7,no'0'
        #print(k)
        counts_all[k] = classmap_list.count(k)
    interest_number_all=counts_all.sum()-counts_all[6]#remove background
    
    for t in ['far','around_inside','whole']:
        locals()[t]=[]
    for class_name in classes_for_analysis:
        locals()['far_'+class_name]=len([n for n in eval(class_name+'_dis') if n < -20])/interest_number_all*100
        locals()['around_'+class_name]=len([n for n in eval(class_name+'_dis') if n >= (-20)and n<=10])/interest_number_all*100
        locals()['inside_'+class_name]=len([n for n in eval(class_name+'_dis') if n>10])/interest_number_all*100
        #merge------------
        locals()['around_inside_'+class_name]=eval('around_'+class_name)+eval('inside_'+class_name)
        locals()['whole_'+class_name]=eval('around_'+class_name)+eval('inside_'+class_name)+eval('far_'+class_name)
        
        for t in ['far','around_inside','whole']:
            locals()[t].append(eval(t+'_'+class_name))
    
    
    #dis=====组间距离和组内距离==================
    areas=['inside_tumor','far_tumor','whole']
    CH_scores=[]
    SDIs=[]
    for area in areas:
        if area=='whole':
            classmap2=classmap
        elif area=='inside_tumor':
            classmap2=classmap*Tumor_mask
        elif area=='far_tumor':
            classmap2=classmap*(1-Tumor_mask)
#        colors=['linen','red','yellow','brown','green','cyan','linen','blue']
#        classmap_rgb = color.label2rgb(classmap,colors=colors,kind='overlay')
#        plt.imshow(classmap_rgb)
#    
        location=[]
        tissue_type=[]
        for i in range(a):
            for j in range(b):
                c=classmap2[i,j]
                if c!=0 and c!=6:
                    location.append([i,j])
                    tissue_type.append(c)
        location=np.float64(np.array(location))
        tissue_type=np.int32(np.array(tissue_type))
        
        CH_score,extra_disp,intra_disp=metrics.calinski_harabasz_score(location, tissue_type)
        CH_scores.append(CH_score)
        #===============================
        SDI_slide=[]
        for i in range(classmap2.shape[0]):
            for j in range(classmap2.shape[1]):
                #print(i,j)
                #i,j=101,72
                if i-5>0 and i+5<a and j-5>0 and j+5<b:
                    submap=classmap2[i-5:i+5,j-5:j+5]
                    submap=submap.tolist()
                    submap = list(itertools.chain(*submap))
                    counts=np.zeros(8)
                    for k in range(1,8):#1~9,no'0'
                        #print(k)
                        counts[k] = submap.count(k)
                        #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6','Mucus==7']
                    #print(counts)
                    All=sum(counts)
                    All_interest=All-counts[6]
                    if All_interest>0.9:
                        counts_SDI=counts[np.array([1,2,3,4,5,7])]
                        probs =counts_SDI / sum(counts_SDI)
                        SDI = stats.entropy(probs, base=2)
                        SDI_slide.append(SDI)
        SDIs.append(np.mean(SDI_slide))
    for i,t in enumerate(['far','around_inside','whole']):
        #print(i,t)
        locals()[t].append(CH_scores[i])
        locals()[t].append(SDIs[i])
    
    classmap2=classmap*Tumor_mask       
    Tumor_index=np.where(classmap2==1)
    Tumor_index=np.array(Tumor_index)
    for class_name in classes_for_analysis:
        locals()[class_name+"_Interaction"]=0
    norm=0
    for index in Tumor_index.transpose().tolist():
        i=index[0]
        j=index[1]
        #print(i,j)
        #i,j=175,60
        
        if i-1>1 and i+2<a and j-1>1 and j+2<b and (i%2==0) and (j%2==0):
            norm=norm+1
            submap=classmap[i-1:i+2,j-1:j+2]
            #print(i,j)
            #print(submap)
            submap=submap.tolist()
            submap = list(itertools.chain(*submap))
            counts=np.zeros(8)
            for k in range(1,8):#1~7,no'0'
                #print(k)
                counts[k] = submap.count(k)
                #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6','Mucus==7']
            for class_name in classes_for_analysis:
                locals()[class_name+"_Interaction"]=eval(class_name+"_Interaction")+counts[class_number_relation_dict[class_name]]
    for class_name in classes_for_analysis:
        locals()[class_name+"_Interaction"]=eval(class_name+"_Interaction")/norm
        #print(eval(class_name+"_Interaction"))
        locals()['around_inside'].append(eval(class_name+"_Interaction"))
                        
    feature_vector=[eval('far'),eval('around_inside'),eval('whole')]
    outputt=[feature_vector,res,LYM_dis,STR_dis]
    #feature_vector=DataFrame({'far':far,'inside':around_inside,'whole':whole})
    return(outputt)
    
#['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]

from concurrent.futures import ProcessPoolExecutor

from tqdm import tqdm
#import concurrent.futures
def run(function, inputs):
    #with concurrent.futures.ThreadPoolExecutor() as executor:
    with ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(function, inputs), total=len(inputs)))
    return results


#global area
#area='whole'
threshold_tumor_size=0
ress=[]
for res in reslist_effect:
    if get_tumor_size(res)>threshold_tumor_size:
        ress.append(res)
    else:
        pass

results=run(calculate_feature2, ress)

All_slides_feture=[]
for result in results:
    if not result==None:
        All_slides_feture.append(result)


feature_name1=['Stroma_proportion','Necrosis_proportion','Lymphocyte_proportion','Normal_hepatic_proportion','Mucus_proportion','CH_score','SDI']
feature_name2=['Stroma_proportion','Necrosis_proportion','Lymphocyte_proportion','Normal_hepatic_proportion','Mucus_proportion','CH_score','SDI','Stroma_interaction','Necrosis_interaction','Lymphocyte_interaction','Normal_hepatic_interaction','Mucus_interaction']

feature_name_all=[]
for i,t in enumerate(['far','infiltrated','whole']):
    if i!=1:
        feature_name_all.append([t+'_'+name for name in feature_name1])
    else:
        feature_name_all.append([t+'_'+name for name in feature_name2])
feature_name_all=list(itertools.chain(*feature_name_all))
All_feature=[]
res_name=[]
All_LYM_dis=[]
All_STR_dis=[]
for fe in All_slides_feture:
    temp=[]
    for i,t in enumerate(['far','infiltrated','whole']):
        temp.extend(fe[0][i])  
    All_feature.append(temp)
    res_name.append(fe[1])
    All_LYM_dis.append(fe[2])
    All_STR_dis.append(fe[3])

import pickle
f=open('distance_ZZ.pkl','wb')
pickle.dump([All_LYM_dis,All_STR_dis],f)
f.close()
'''
feature = feature_name_all
name=feature
test=pd.DataFrame(columns=name,data=All_feature)#
print(test)
test.index =res_name#[res.split('.')[0] for res in reslist]
test.to_csv('/data/backup/Yuni/CRC_Liver/12_feature_calculation/01/all_feature_new.csv',encoding='gbk')
'''


All_LYM_dis=list(itertools.chain.from_iterable(All_LYM_dis))
All_STR_dis=list(itertools.chain.from_iterable(All_STR_dis))
#------------画分布直方图---------------

plt.figure()
plt.title("LYM distance to tumor")
sns.set_style('darkgrid')
Blue_dia_all1= list(itertools.chain.from_iterable(All_LYM_dis))
sns.distplot(Blue_dia_all1)#,hist_kws={'color':'green'}





#-----------------高斯混合分布-------------------------
def Normal(x,mu,sigma):#一元正态分布概率密度函数
    
    return np.exp(-(x-mu)**2/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)

def GMM2(data,n,Mu,SigmaSquare,Alpha):
    N=data.shape[0]#观测集大小
    #n:迭代次数   
    np.random.seed(1)
    
    i=0#迭代次数
    
    while(i<=n):#用EM算法迭代求参数估计
        i+=1
        #Expectation--------
        if SigmaSquare[0][0]<=0.001:
            SigmaSquare[0][0]=SigmaSquare[0][0]+0.1
        if SigmaSquare[0][1]<=0.001:
            SigmaSquare[0][1]=SigmaSquare[0][1]+0.1
        gauss1=Normal(data,Mu[0][0],np.sqrt(SigmaSquare[0][0]))#第一个模型
        gauss2=Normal(data,Mu[0][1],np.sqrt(SigmaSquare[0][1]))#第二个模型
        
        Gamma1=Alpha[0][0]*gauss1
        Gamma2=Alpha[0][1]*gauss2
        M=Gamma1+Gamma2+0.00001

    
        #Gamma=np.concatenate((Gamma1/m,Gamma2/m),axis=1) 元素(j,k)为第j个样本来自第k个模型的概率，聚类时用来判别样本分类
    
        #Maximization--------
        #更新SigmaSquare
        SigmaSquare[0][0]=np.dot((Gamma1/M).T,(data-Mu[0][0])**2)/np.sum(Gamma1/M)
        SigmaSquare[0][1]=np.dot((Gamma2/M).T,(data-Mu[0][1])**2)/np.sum(Gamma2/M)
    
        #更新Mu  
        Mu[0][0]=np.dot((Gamma1/M).T,data)/np.sum(Gamma1/M)
        Mu[0][1]=np.dot((Gamma2/M).T,data)/np.sum(Gamma2/M)
    
        #更新Alpha
        Alpha[0][0]=np.sum(Gamma1/M)/N
        Alpha[0][1]=np.sum(Gamma2/M)/N
        
        if(i%1==0):
            print ("第",i,"次迭代:")
            print ("Mu:",Mu)
            print ("Sigma:",np.sqrt(SigmaSquare))
            print ("Alpha",Alpha)
    return(Mu,SigmaSquare,Alpha)


def GMM3(data,n,Mu,SigmaSquare,Alpha):

    N=data.shape[0]#观测集大小
    #n:迭代次数
    np.random.seed(1)

    i=0#迭代次数
    
    while(i<=n):#用EM算法迭代求参数估计
        
        i+=1
        
        #Expectation
        if SigmaSquare[0][0]<=0.001:
            SigmaSquare[0][0]=SigmaSquare[0][0]+0.03
        if SigmaSquare[0][1]<=0.001:
            SigmaSquare[0][1]=SigmaSquare[0][1]+0.03
        if SigmaSquare[0][2]<=0.001:
            SigmaSquare[0][2]=SigmaSquare[0][2]+0.03
        gauss1=Normal(data,Mu[0][0],np.sqrt(SigmaSquare[0][0]))#第一个模型
        gauss2=Normal(data,Mu[0][1],np.sqrt(SigmaSquare[0][1]))#第二个模型
        gauss3=Normal(data,Mu[0][2],np.sqrt(SigmaSquare[0][2]))#第三个模型
        
        Gamma1=Alpha[0][0]*gauss1
        Gamma2=Alpha[0][1]*gauss2
        Gamma3=Alpha[0][2]*gauss3
    
        M=Gamma1+Gamma2+Gamma3+0.00001
    
        #Gamma=np.concatenate((Gamma1/m,Gamma2/m),axis=1) 元素(j,k)为第j个样本来自第k个模型的概率，聚类时用来判别样本分类
    
        #Maximization
        
        #更新SigmaSquare
        
        SigmaSquare[0][0]=np.dot((Gamma1/M).T,(data-Mu[0][0])**2)/np.sum(Gamma1/M)
        
        SigmaSquare[0][1]=np.dot((Gamma2/M).T,(data-Mu[0][1])**2)/np.sum(Gamma2/M)
        SigmaSquare[0][2]=np.dot((Gamma3/M).T,(data-Mu[0][2])**2)/np.sum(Gamma3/M)
    
        #更新Mu       
    
        Mu[0][0]=np.dot((Gamma1/M).T,data)/np.sum(Gamma1/M)
        Mu[0][1]=np.dot((Gamma2/M).T,data)/np.sum(Gamma2/M)
        Mu[0][2]=np.dot((Gamma3/M).T,data)/np.sum(Gamma3/M)
    
        #更新Alpha
    
        Alpha[0][0]=np.sum(Gamma1/M)/N  
        Alpha[0][1]=np.sum(Gamma2/M)/N
        Alpha[0][2]=np.sum(Gamma3/M)/N
        
        if(i%10==0):
            print ("第",i,"次迭代:")
            print ("Mu:",Mu)
            print ("Sigma:",np.sqrt(SigmaSquare))
            print ("Alpha",Alpha)
    return(Mu,SigmaSquare,Alpha)


def gaussian(sigma, x, u):
	y = np.exp(-(x - u) ** 2 / (2 * sigma ** 2)) / (sigma * math.sqrt(2 * math.pi))
	return y

#蓝色细胞核大小-----------
sns.set_style('darkgrid')
Blue_dia_all1=list(itertools.chain.from_iterable(All_LYM_dis))
Blue_dia_all_array=np.array(Blue_dia_all1)
Blue_dia_all_array.shape=Blue_dia_all_array.shape[0],1
data=Blue_dia_all_array#合并身高数据，N行1列
#data=data[np.where(data < 25)]
#data=data[np.where(data > 5)]


#GMM2
Mu=np.array([[-50,0.5]])#平均值向量
SigmaSquare=np.array([[2,8]]) #模型迭代用Sigma平方
Alpha=np.array([[0.5,0.5]])  
n=300
Mu,SigmaSquare,Alpha=GMM2(data,n,Mu,SigmaSquare,Alpha)
x = np.linspace(data.min(), data.max(), 1000)
plt.figure()
plt.title("LYM distance to tumor_GMM2")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1)
a=Mu[0][1]
b=gaussian(np.sqrt(SigmaSquare[0][1]), Mu[0][1], Mu[0][1])*Alpha[0][1]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.show()



#GMM3
Mu=np.array([[-75,-50,20]])#平均值向量
SigmaSquare=np.array([[2.0,6.5,8.0]]) #模型迭代用Sigma平方
Alpha=np.array([[0.33,0.33,0.34]])#随机初始化各模型比重系数（大于等于0，且和为1）
n=60
Mu,SigmaSquare,Alpha=GMM3(data,n,Mu,SigmaSquare,Alpha)#data,iterate number
x = np.linspace(data.min(), data.max(), 1000)
plt.figure()
plt.title("LYM distance to tumor_GMM3")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][2]), x, Mu[0][2])*Alpha[0][2], "b-", linewidth=1)
a=Mu[0][1]
b=gaussian(np.sqrt(SigmaSquare[0][1]), Mu[0][1], Mu[0][1])*Alpha[0][1]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.show()



All_STR_dis


sns.set_style('darkgrid')
Blue_dia_all1= list(itertools.chain.from_iterable(All_STR_dis))
Blue_dia_all_array=np.array(Blue_dia_all1)
Blue_dia_all_array.shape=Blue_dia_all_array.shape[0],1
data=Blue_dia_all_array#合并身高数据，N行1列
#data=data[np.where(data < 25)]
#data=data[np.where(data > 5)]


#GMM2
Mu=np.array([[-50,0.5]])#平均值向量
SigmaSquare=np.array([[2,8]]) #模型迭代用Sigma平方
Alpha=np.array([[0.5,0.5]])  
n=300
Mu,SigmaSquare,Alpha=GMM2(data,n,Mu,SigmaSquare,Alpha)
x = np.linspace(data.min(), data.max(), 1000)
plt.figure()
plt.title("STR distance to tumor_GMM2")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1)
a=Mu[0][1]
b=gaussian(np.sqrt(SigmaSquare[0][1]), Mu[0][1], Mu[0][1])*Alpha[0][1]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.show()



#GMM3
Mu=np.array([[-75,-50,20]])#平均值向量
SigmaSquare=np.array([[2.0,6.5,8.0]]) #模型迭代用Sigma平方
Alpha=np.array([[0.33,0.33,0.34]])#随机初始化各模型比重系数（大于等于0，且和为1）
n=100
Mu,SigmaSquare,Alpha=GMM3(data,n,Mu,SigmaSquare,Alpha)#data,iterate number
x = np.linspace(data.min(), data.max(), 1000)
plt.figure()
plt.title("STR distance to tumor_GMM3")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][2]), x, Mu[0][2])*Alpha[0][2], "b-", linewidth=1)
a=Mu[0][1]
b=gaussian(np.sqrt(SigmaSquare[0][1]), Mu[0][1], Mu[0][1])*Alpha[0][1]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.show()

'''


















'''
ALL_Tum_score=[]
Entropy_slides=[]
ALL_Tum_intereaction_score=[]
Disp_slides=[]
ALL_LYM_scores=[]
ALL_STR_scores=[]
result_success=[]
for res in reslist_chosen:
#    try:
    #res='754206_10B'
    print('processing:',res)
        #res=reslist1_new[0]
        
    f=open(f'{path_res}{res}.pkl','rb')
    #name,predicted_classes,label,classmap,x,y,classes=pickle.load(f)
    name,predicted_classes,label=pickle.load(f)
    #Tumscore_all.append(Tum_score)
    #classmap2=np.array(classmap)
    #colors=['linen','red','yellow','brown','green','cyan','linen','blue']
    #colors =[(0.41,0.41,0.41), (1,1,1), (0.33,0.25,0.27), (0.76470588, 0.39215686, 0.77254902),(0.98823529, 0.42352941, 0.52156863), (0.80392157, 0.36078431, 0.36078431), (1, 0.63921569, 0.2627451), (0.2745098 , 0.50980392, 0.70588235), (0.10980392, 0.6745098 , 0.47058824)]#,
#    test=np.array([[0,1,2,3,4,5,6,7],[0,1,2,3,4,5,6,7],[0,1,2,3,4,5,6,7]])
#    classmap_test = color.label2rgb(test,colors=colors,kind='overlay')
#    plt.imshow(classmap_test)
    #classmap_rgb = color.label2rgb(classmap2,colors=colors,kind='overlay')
    f.close()
#    plt.style.use('classic')
#    plt.figure()
#    plt.imshow(classmap_rgb)
    slideImage = openslide.open_slide(f'{path_svs}{res}.svs')
    dim=slideImage.level_dimensions[0]#1
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
    
    a,b=classmap.shape
    
#    dst=transform.resize(classmap, (a*224, b*224))
#    plt.imshow(dst)
    #寻找肿瘤区域，并计算STR和LYM到肿瘤区域的距离分布
    Tumor_mask=(classmap==1)
    #plt.imshow(Tumor_mask)
    Tumor_mask=Tumor_mask.astype(np.uint8)*255
    kernel = morphology.disk(5)
    #kernel2 = np.ones((5,5),np.uint8) 
    Tumor_mask = cv2.morphologyEx(Tumor_mask, cv2.MORPH_CLOSE, kernel)
    Tumor_mask = cv2.morphologyEx(Tumor_mask, cv2.MORPH_CLOSE, kernel)
    #Tumor_mask =morphology.closing(Tumor_mask, selem=None, out=None)
    #Tumor_mask=morphology.remove_small_holes(Tumor_mask, 1000)
    Tumor_mask=ndimage.binary_fill_holes(Tumor_mask)
    #Tumor_mask=(Tumor_mask==255)
    Tumor_mask=morphology.remove_small_objects(Tumor_mask,50)
    #plt.imshow(Tumor_mask)
    Tumor_mask=Tumor_mask.astype(np.uint8)*255

    contours, hierarchy = cv2.findContours(Tumor_mask, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    
    STR_dis=[]
    LYM_dis=[]
    
#        approx = cv2.approxPolyDP(cnt, 3, True)
#        # 3.画出多边形
#        Tumor_mask = cv2.cvtColor(Tumor_mask, cv2.COLOR_GRAY2BGR)
#        cv2.polylines(Tumor_mask, [approx], True, (0, 255, 0), 2)
    for i in range(a):
        for j in range(b):
            c=classmap[i,j]
            if c==2 or c==4:
                temp=[]
                for cnt in contours:
                    temp.append(cv2.pointPolygonTest(cnt, (i, j), True))
                dis=max(temp)
                if c==2:
                    STR_dis.append(dis)
                else:
                    LYM_dis.append(dis)
    
    
    classmap_list=classmap.tolist()
    classmap_list = list(itertools.chain(*classmap_list))
    counts_all=np.zeros(8)
    for k in range(1,8):#1~9,no'0'
        #print(k)
        counts_all[k] = classmap_list.count(k)
    interest_number_all=counts_all[1]#+counts_all[2]+counts_all[3]+counts_all[4]+counts_all[7]
    
    far_LYM=len([n for n in LYM_dis if n < -20])/interest_number_all*100
    around_LYM=len([n for n in LYM_dis if n >= (-20)and n<=0])/interest_number_all*100
    inside_LYM=len([n for n in LYM_dis if n>0])/interest_number_all*100
    #near_LYM=around_LYM+inside_LYM??
    #ALL_LYM_scores.append([far_LYM,around_LYM,inside_LYM])
    far_STR=len([n for n in STR_dis if n < -20])/interest_number_all*100
    around_STR=len([n for n in STR_dis if n >= (-20)and n<=0])/interest_number_all*100
    inside_STR=len([n for n in STR_dis if n>0])/interest_number_all*100
    #ALL_STR_scores.append([far_STR,around_STR,inside_STR])
    
    #sns.set_style('darkgrid')
    #sns.set_context(rc={'figure.figsize': (8, 5) } )
    #plt.hist(LYM_dis)
    #plt.hist(STR_dis)
    #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]
    #dis=====组间距离和组内距离==================
    location=[]
    tissue_type=[]
    for i in range(a):
        for j in range(b):
            c=classmap[i,j]
            if c!=0 and c!=6:
                location.append([i,j])
                tissue_type.append(c)
    location=np.float64(np.array(location))
    tissue_type=np.int32(np.array(tissue_type))
    
    CH_score,extra_disp,intra_disp=metrics.calinski_harabasz_score(location, tissue_type)
    #Disp_slides.append([extra_disp,intra_disp,CH_score])
    #===============================
    
    #==========比例以及和肿瘤相邻的比例===========================
    
    #plt.imshow(Tumor_mask)
    #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]
    Tumor_index=np.where(classmap==1)
    
    random.seed(10)
    sample_N=len(Tumor_index[1])
    num=list(range(0,sample_N))
    random_num = random.sample(num, sample_N)
    Sample_num=int(sample_N*0.5)
    Tumor_index=np.array(Tumor_index)  
    random_index=Tumor_index[:,random_num[1:Sample_num]]
    #STR_mask=(classmap==2)
    #plt.imshow(STR_mask)
    STR_Tum_slide=[]
    LYM_Tum_slide=[]
    NORM_Tum_slide=[]
    DEB_Tum_slide=[]
    
    STR_Tum_Interaction_slide=[]
    LYM_Tum_Interaction_slide=[]

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
            counts=np.zeros(7)
            for k in range(1,7):#1~9,no'0'
                #print(k)
                counts[k] = submap.count(k)#['BACK==0','ADI','BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']
                #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK']
            #print(counts)
            All=sum(counts)
            All_interest=All-counts[6]
            if counts[1]/All>=0.2 and counts[1]/All<0.85:#
                STR_Tum=counts[2]/All_interest
                STR_Tum_slide.append(STR_Tum)
                #print(STR_Tum_slide)
                LYM_Tum=counts[4]/All_interest
                LYM_Tum_slide.append(LYM_Tum)
                NORM_Tum=counts[5]/All_interest
                NORM_Tum_slide.append(NORM_Tum)
                DEB_Tum=counts[3]/All_interest
                DEB_Tum_slide.append(DEB_Tum)

                STR_Tum_Interaction=0
                LYM_Tum_Interaction=0
                
                submap=classmap[i-5:i+5,j-5:j+5]
                Tumor_index2=np.where(submap==1)
                Tumor_index2=np.array(Tumor_index2) 
                for index2 in Tumor_index2.transpose().tolist():
                    m=index2[0]
                    n=index2[1]
                    #print(m,n)
                    if m-1>=0 and m+1<11 and n-1>=0 and n+1<11 and (m%3==0):
                        submap2=submap[m-1:m+2,n-1:n+2]
                        #print("-----index:",m,n)
                        #print(submap)
                        submap2=submap2.tolist()
                        submap2 = list(itertools.chain(*submap2))
                        counts2=np.zeros(7)
                        for k in range(1,7):#1~9,no'0'
                            #print(k)
                            counts2[k] = submap2.count(k)#['BACK==0','ADI','BACK', 'DEB', 'LYM', 'MUC', 'MUS', 'NORM', 'STR', 'TUM']
                            #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK']
                        del submap2
                        #print("counts:",counts2)
                        #print("STR:",counts2[2])
                        STR_Tum_Interaction=STR_Tum_Interaction+counts2[2]
                        LYM_Tum_Interaction=LYM_Tum_Interaction+counts2[4]

                #print("STR_Tum_Interaction:",STR_Tum_Interaction)    
                #print("LYM_Tum_Interaction:",LYM_Tum_Interaction)        
                
                STR_Tum_Interaction_slide.append(STR_Tum_Interaction/All_interest)
                LYM_Tum_Interaction_slide.append(LYM_Tum_Interaction/All_interest)
                
    try:
        All_proportion_slide=[STR_Tum_slide,LYM_Tum_slide,DEB_Tum_slide,NORM_Tum_slide]
        Tum_score=[np.percentile(pro,90) for pro in All_proportion_slide]
        ALL_Tum_score.append(Tum_score)
        Tum_score=np.array(Tum_score)
    #=====================================
    
    #plt.pause(100)
    
    
        All_Interaction_proportion_slide=[STR_Tum_Interaction_slide,LYM_Tum_Interaction_slide]
    
        Tum_intereaction_score=[np.percentile(pro,90) for pro in All_Interaction_proportion_slide]
    except:
        print(f'{res} failed!')
        continue
    ALL_Tum_intereaction_score.append(Tum_intereaction_score)
    #计算概率分布==========计算熵=========
    probs =Tum_score / sum(Tum_score)
    S = stats.entropy(probs, base=2)
    Entropy_slides.append(S)
    result_success.append(res.split('.')[0])
    Disp_slides.append([extra_disp,intra_disp,CH_score])
    ALL_LYM_scores.append([far_LYM,around_LYM,inside_LYM])
    ALL_STR_scores.append([far_STR,around_STR,inside_STR])
#    except:
#        print(res,"failed!!!")
#        continue
#====boxplot==============================
#    all_data = All_Interaction_proportion_slide
#    labels = ['STR_TUM', 'LYM_TUM']
#    
#    bplot = plt.boxplot(all_data, patch_artist=True, labels=labels)  # 设置箱型图可填充
#    plt.title('Tumor Interaction')
#    
#    colors = ['pink', 'lightblue']
#    for patch, color1 in zip(bplot['boxes'], colors):
#        patch.set_facecolor(color1)  # 为不同的箱型图填充不同的颜色
#    
#    #plt.yaxis.grid(True)
#    plt.xlabel('Interaction Types')
#    plt.ylabel('Interation Level')
#    plt.show()

#============radar graph==========================
#
#    plt.rcParams['font.sans-serif'] = 'Microsoft YaHei'
#    plt.rcParams['axes.unicode_minus'] = False
#    
#
#    plt.style.use('ggplot')
#    
#
#    values = Tum_score
#    feature = ['STR','LYM','DEB','NORM']
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
#    ax.set_ylim(0,0.8)
#
#    plt.title('Cell Proportion')
#
#    ax.grid(True)
#    #plt.savefig(f'{outPath}_interaction.jpg')
#
#    plt.show()
        

#np.savetxt('Tumscore_all.csv', Tumscore_all, delimiter = ',')
  
if not os.path.exists(f'{root_path}03_clinical_info/'):
       os.makedirs(f'{root_path}03_clinical_info/')  

#
feature = ['STR','LYM','DEB','NORM']
name=feature
test=pd.DataFrame(columns=name,data=ALL_Tum_score)#
print(test)
test.index =result_success
test.to_csv(f'{root_path}03_clinical_info/tum_interection.csv',encoding='gbk')

#
feature = ['Entropy of slide']
name=feature
test=pd.DataFrame(columns=name,data=Entropy_slides)#
print(test)
test.index =result_success
test.to_csv(f'{root_path}03_clinical_info/tum_Entropy.csv',encoding='gbk')
#
##
feature = ['STR','LYM']
name=feature
test=pd.DataFrame(columns=name,data=ALL_Tum_intereaction_score)#
print(test)
test.index =result_success
test.to_csv(f'{root_path}03_clinical_info/ALL_Tum_intereaction_score.csv',encoding='gbk')
#
##
feature = ['Extra_disp','Intra_disp','CH_score']
name=feature
test=pd.DataFrame(columns=name,data=Disp_slides)#
print(test)
test.index =result_success
test.to_csv(f'{root_path}03_clinical_info/Extra_disp_score.csv',encoding='gbk')

##
feature = ['far_LYM','around_LYM','inside_LYM']
name=feature
test=pd.DataFrame(columns=name,data=ALL_LYM_scores)#
print(test)
test.index =result_success
test.to_csv(f'{root_path}03_clinical_info/distance_LYM_score.csv',encoding='gbk')

#
feature = ['far_STR','around_STR','inside_STR']
name=feature
test=pd.DataFrame(columns=name,data=ALL_STR_scores)#
print(test)
test.index =result_success#[res.split('.')[0] for res in reslist]
test.to_csv(f'{root_path}03_clinical_info/distance_STR_score.csv',encoding='gbk')


'''
