# -*- coding: utf-8 -*-
import os
import pickle
import numpy as np
from skimage import color
import matplotlib.pyplot as plt
import random
import pandas as pd
import itertools
from skimage import morphology
from skimage import transform#,data
import cv2
from scipy import stats
from skimage import measure
from sklearn import metrics
import seaborn as sns
from scipy import ndimage
import math
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
Blue_dia_all1= list(itertools.chain.from_iterable(All_LYM_dis))
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
