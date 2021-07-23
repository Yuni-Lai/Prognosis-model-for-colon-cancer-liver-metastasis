#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import xlrd
import os
import pickle
import numpy as np
from skimage import color
import matplotlib.pyplot as plt
import random
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

from pandas.core.frame import DataFrame

#clinical=pd.read_excel("/data/backup/Yuni/CRC_Liver/06_clinical_info/20201031_ZZ_CRLM_final_Eng.xlsx",sheet_name='new_145cases')#399cases
###
#f=open('/data/backup/Yuni/CRC_Liver/06_clinical_info/new145_liver_survival.pkl','wb')
#pickle.dump([clinical],f)
#f.close()


f=open('/data/backup/Yuni/CRC_Liver/06_clinical_info/new145_liver_survival.pkl','rb')
clinical=pickle.load(f)
f.close()
clinical=clinical[0]
liver_name=[str(name) for name in clinical.iloc[:,1]]

#--------setting---path----------------
#01
#root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_single_85/'
#path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/01_single_85/01_classified_first/'
#path_svs='/data/backup/Yuni/CRC_Liver/00_SVS_images/01_SVS_image_85/'

#02
#root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/02_single_285/'
#path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/02_single_285/01_classified_second/'
#path_svs='/data/backup/Yuni/CRC_Liver/00_SVS_images/02_liver_SVS_285/'

#03
#global path_res
#global path_svs
#global root_path
#root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/03_mult_265/'
#path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/03_mult_265/01_classified_third/'
#path_svs='/data/backup/Yuni/CRC_Liver/00_SVS_images/03_multi_regions_SVS_265/'

#04
global path_res
global path_svs
global root_path
root_path='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/04_mult_141/'
path_res='/data/backup/Yuni/CRC_Liver/10_Finetuned_Oct/04_mult_141/01_classified/'
path_svs='/data/backup2/Yuni/CRC_Liver/00_SVS_images/04_multi_regions/CRC_LM_10/'

liver_all=pd.read_excel("/data/backup/Yuni/CRC_Liver/06_clinical_info/image_list＿04_info.xlsx",sheet_name='All')#399cases
conlon=pd.read_excel("/data/backup/Yuni/CRC_Liver/06_clinical_info/image_list＿04_info.xlsx",sheet_name='Colon')
colon_number=[str(name).replace(' ','_') for name in conlon.iloc[:,2]]
liver_chosen=liver_all[liver_all["For_single_region_analysis"]=='Yes']#144
liver_chosen_number=[str(name).replace(' ','_') for name in liver_chosen.iloc[:,2]]
len(liver_chosen_number)


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

#reslist = os.listdir(path_res)
reslist = get_image_names(path_res)

reslist_new_number=[res.split('.')[0] for res in reslist]

reslist_effect=[]
print("chosen for sing region analysis:------")
for res in reslist_new_number:
    if res.split('/')[1] in [str(name).split('.')[0] for name in liver_chosen_number]:
        reslist_effect.append(res)
        print(res)
print("-------------------------")
#Tumscore_all=[]
#classmap_rgb_all=[]

#import skimage

res_patients=[res.split('_')[0] for res in reslist]
patients=dict(Counter(res_patients))

#['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]
def get_tumor_size(res):
    #res=reslist[0].split('.')[0]
    #print('processing:',res)
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



def calculate_feature2(inputt):#according to tumor region
#    try:
    #res=reslist[20].split('.')[0]
    #res='2020-06-29/746472_1G,1H'
    patient=inputt[0]
    res=inputt[1]
    print('processing:',res)
        #res=reslist1_new[0]
    try:
        f=open(f'{path_res}{res}.pkl','rb')
        #name,predicted_classes,label,classmap,x,y,classes=pickle.load(f)
        name,predicted_classes,label=pickle.load(f)
        f.close()
    except:
        return None
    if res =='2020-07-10/592136_1A':
        res = '2020-07-10/592136_1A,1B.svs'

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
        #a,b=n.split('_')
        a,b=n[0],n[1]
        a,b=int(a),int(b)
        classmap[b,a]=predicted_classes[i]+1#['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK']=1~9; 0=white back ground
        x[i],y[i]=b,a
    if index_max==0:      
        classmap=classmap[0:a2,:]
    else:
        classmap=classmap[:,0:a1]
    
    a,b=classmap.shape
    #plt.imshow(classmap)
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
        locals()['far_'+class_name]=len([n for n in eval(class_name+'_dis') if n < -10])/interest_number_all*100
        locals()['around_'+class_name]=len([n for n in eval(class_name+'_dis') if n >= (-10)and n<=0])/interest_number_all*100
        locals()['inside_'+class_name]=len([n for n in eval(class_name+'_dis') if n>0])/interest_number_all*100
        #merge------------
        locals()['around_inside_'+class_name]=eval('around_'+class_name)+eval('inside_'+class_name)
        locals()['whole_'+class_name]=eval('around_'+class_name)+eval('inside_'+class_name)+eval('far_'+class_name)
        
        for t in ['far','around_inside','whole']:
            locals()[t].append(eval(t+'_'+class_name))
    
    
    
    
    #sns.set_style('darkgrid')
    #sns.set_context(rc={'figure.figsize': (8, 5) } )
    #plt.hist(LYM_dis)
    #plt.hist(STR_dis)
    #['BACK==0','Tumor==1','Stroma==2','Necrosis==3','Lymphocyte==4','Normal_hepatic==5','BACK==6'，Mucus==7]
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
    outputt=[patient,feature_vector,res]
    #feature_vector=DataFrame({'far':far,'inside':around_inside,'whole':whole})
    return(outputt)


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
thread_num=1
thread_tumor_size=8000

All_slides_feture=[]
reslist_chosen=list(patients.keys())
ress=[]
for i,p in enumerate(patients):
    print(p,'----------')
    feature_one_patient=[]
    ress_temp=[]
    
    for res in reslist_new_number:
        if res.split('_')[0]==p:
            if res.split('/')[1]+'.svs' not in colon_number:
                if get_tumor_size(res)>thread_tumor_size:
#                    feature_one_slide=calculate_feature(res,'whole')
#                    feature_one_patient.append(feature_one_slide)
                    ress_temp.append([p,res])
                    #print(ress_temp)
    if len(ress_temp)>thread_num:
        ress.extend(ress_temp)
    else:
        pass

len(ress)
results=run(calculate_feature2, ress)


pp=[]#patient that calculated
for result in results:
    pp.append(result[0])

pp=dict(Counter(pp))
pp=list(pp.keys())
All_slides_feture=[]            
for p in pp:
    feature_one_patient=[]
    for result in results:
        if result[0]==p:
            feature_one_patient.append(result[1])
    All_slides_feture.append(feature_one_patient)

#f=open('All_slides_feture.pkl','wb')
#pickle.dump(All_slides_feture,f)
#f.close()
#
#
#f=open('results.pkl','wb')
#pickle.dump(results,f)
#f.close()


#area='inside_tumor'
#inside_tumor_feture=[]
#reslist_chosen=list(patients.keys())
#for i,p in enumerate(patients):
#    print(p,'----------')
#    feature_one_patient=[]
#    ress=[]
#    for res in reslist_new_number:
#        if res.split('_')[0]==p:
#            if res.split('/')[1]+'.svs' not in colon_number:
#                if get_tumor_size(res)>300:
#                    #feature_one_slide=calculate_feature(res,'inside_tumor')
#                    #feature_one_patient.append(feature_one_slide)
#                    ress.append(res)
#    if len(ress)>3:
#        results=run(calculate_feature, ress)
#        for result in results:
#            feature_one_patient.append(result)
#        inside_tumor_feture.append(feature_one_patient)
#    else:
#        pass
    


def test_corr(df,plot,t):
    df_norm=(df-eval(t+'_Min'))/(eval(t+'_Max')-eval(t+'_Min'))

    #df_norm=(df - eval(t+'_Mean')) / (eval(t+'_Std'))
    df_norm=df_norm.T

    
    dfData_corr = df_norm.corr(method='spearman')#'pearson',spearman
    #dfData_corr =df.corr()
    mean_Cor=np.mean((dfData_corr-np.identity(dfData_corr.shape[0])).mean())
    mean_Cor=round(mean_Cor, 3)
    if plot==1:
        plt.subplots(figsize=(9, 9)) # 设置画面大小
        plt.title(f'{t}, mean of Cor:{mean_Cor}')
        sns.heatmap(dfData_corr, annot=True, vmax=1, square=True, cmap="Blues")
    #plt.savefig('./BluesStateRelation.png')
        plt.show()
        plt.subplots(figsize=(15, 9)) # 设置画面大小
        plt.title(f'feature of {t}')
        sns.heatmap(df_norm, annot=True, vmax=1, square=False, cmap="Greens")
        plt.show()
    else:
        pass
    return(mean_Cor)



types=['far','around_inside','whole']#
for i,t in enumerate(types):
    F=[]
    for feature in All_slides_feture:
        #feature=All_slides_feture[0]
        if len(feature)>thread_num:
            #t='around_inside'
            #i=1    
            for f in feature:
                F.append(list(f[i][0:7]))
    df=DataFrame(np.array(F))
    locals()[t+'_Mean']=df.mean()
    locals()[t+'_Std']=df.std()
    locals()[t+'_Min']=df.min()
    locals()[t+'_Max']=df.max()
            
cor_all_type=[]           
types=['far','around_inside','whole']#
for i,t in enumerate(types):
    cor_mean_all_slides=[]
    for k,feature in enumerate(All_slides_feture):
        #feature=All_slides_feture[0]
        if len(feature)>thread_num:
            #t='around_inside'
            #i=1
            F=[]        
            for f in feature:
                #f=feature[0]
                F.append(list(f[i][0:7]))#
            df=DataFrame(np.array(F))
            if k==0:
                flag=1
            else:
                flag=0
            cor_mean=test_corr(df,flag,t)
            cor_mean_all_slides.append(cor_mean)
    cor_all_type.append(cor_mean_all_slides)

plt.subplots(figsize=(9, 9))
plt.title('Boxplot of correlation')
df1=DataFrame(np.array(cor_all_type).T)
df1.columns = types
ax = sns.boxplot(data=df1,dodge=False).set(ylabel='Mean of correlation')#x="types", hue="types"
df1.mean()


area='inside_tumor'
all_res=[]#patient that calculated
for result in results[0:4]:
    all_res.append(result[2])
    f=open(f'{path_res}{result[2]}.pkl','rb')
    #name,predicted_classes,label,classmap,x,y,classes=pickle.load(f)
    name,predicted_classes,label=pickle.load(f)
    f.close()
    slideImage = openslide.open_slide(f'{path_svs}{result[2]}.svs')
    dim=slideImage.level_dimensions[1]#1
    a1=int(dim[0]/224)*2+1
    a2=int(dim[1]/224)*2+1
    maxxx=max(a1,a2)
    index_max=np.argmax(dim)
    classmap=np.zeros([maxxx,maxxx])
    x=np.zeros(len(name))
    y=np.zeros(len(name))
    for i,n in enumerate(name):#'0_0'
        #a,b=n.split('_')
        a,b=n[0],n[1]
        a,b=int(a),int(b)
        classmap[b,a]=predicted_classes[i]+1#['Tumor','Stroma','Necrosis','Lymphocyte','Normal_hepatic','BACK']=1~9; 0=white back ground
        x[i],y[i]=b,a
    if index_max==0:      
        classmap=classmap[0:a2,:]
    else:
        classmap=classmap[:,0:a1]
    
    a,b=classmap.shape
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
    if area=='whole':
        classmap2=classmap
    elif area=='inside_tumor':
        classmap2=classmap*Tumor_mask
    elif area=='far_tumor':
        classmap2=classmap*(1-Tumor_mask)
    contours, hierarchy = cv2.findContours(Tumor_mask, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)#Tumo
    
    plt.figure()
    colors=['linen','red','yellow','brown','green','cyan','linen','blue']
    classmap_rgb = color.label2rgb(classmap2,colors=colors,kind='overlay')
    plt.imshow(classmap_rgb)

#
#cor_mean_inside_tumor_slides=[]
#for feature in inside_tumor_feture:
#    if len(feature)>3:
#        df=DataFrame(np.array(feature).T)
#        cor_mean=test_corr(df,0)
#        cor_mean_inside_tumor_slides.append(cor_mean)       
#np.mean(cor_mean_inside_tumor_slides)





































#np.savetxt('Tumscore_all.csv', Tumscore_all, delimiter = ',')
  
#if not os.path.exists(f'{root_path}03_clinical_info/'):
#       os.makedirs(f'{root_path}03_clinical_info/')  
#
##
#feature = ['STR','LYM','DEB','NORM']
#name=feature
#test=pd.DataFrame(columns=name,data=ALL_Tum_score)#
#print(test)
#test.index =result_success
#test.to_csv(f'{root_path}03_clinical_info/tum_interection.csv',encoding='gbk')
#
##
#feature = ['Entropy of slide']
#name=feature
#test=pd.DataFrame(columns=name,data=Entropy_slides)#
#print(test)
#test.index =result_success
#test.to_csv(f'{root_path}03_clinical_info/tum_Entropy.csv',encoding='gbk')
##
###
#feature = ['STR','LYM']
#name=feature
#test=pd.DataFrame(columns=name,data=ALL_Tum_intereaction_score)#
#print(test)
#test.index =result_success
#test.to_csv(f'{root_path}03_clinical_info/ALL_Tum_intereaction_score.csv',encoding='gbk')
##
###
#feature = ['Extra_disp','Intra_disp','CH_score']
#name=feature
#test=pd.DataFrame(columns=name,data=Disp_slides)#
#print(test)
#test.index =result_success
#test.to_csv(f'{root_path}03_clinical_info/Extra_disp_score.csv',encoding='gbk')
#
###
#feature = ['far_LYM','around_LYM','inside_LYM']
#name=feature
#test=pd.DataFrame(columns=name,data=ALL_LYM_scores)#
#print(test)
#test.index =result_success
#test.to_csv(f'{root_path}03_clinical_info/distance_LYM_score.csv',encoding='gbk')
#
##
#feature = ['far_STR','around_STR','inside_STR']
#name=feature
#test=pd.DataFrame(columns=name,data=ALL_STR_scores)#
#print(test)
#test.index =result_success#[res.split('.')[0] for res in reslist]
#test.to_csv(f'{root_path}03_clinical_info/distance_STR_score.csv',encoding='gbk')
#


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
