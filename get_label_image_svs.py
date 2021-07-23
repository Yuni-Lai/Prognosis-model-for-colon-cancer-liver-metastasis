#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 17:25:40 2020

@author: Yuni
"""

import os
import time

import matplotlib.pyplot as plt
import openslide

#from deepslide_code.utils import get_image_names

start = time.time()
# input_folder = 'C:/Users/Administrator/Desktop/UC/1/'
# output_folder = 'C:/Users/Administrator/Desktop/UC/1/'

input_folder = '/data/backup2/Yuni/CRC_Liver/00_SVS_images/05_beij_511/2020-10-15/'
output_folder = '/data/backup2/Yuni/CRC_Liver/00_SVS_images/05_beij_511/label/'

image_names = os.listdir(input_folder)

for image_name in image_names:
    image_path = os.path.join(input_folder, image_name)
    slide = openslide.open_slide(image_path)
    print("Start: ", image_name)
    outPath = (output_folder + image_name[0:-4] + '.jpg')

    slide_label = slide.associated_images.items()
    a = slide.associated_images['label']
    plt.axis('off')
    plt.imshow(a)
    # fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))
    # for i, (k, v) in enumerate(slide.associated_images.items()):
    #     # print("i=",i)
    #     # print("v=",v)
    #     print("k=",k)
    #     axs[i].imshow(v)
    #     axs[i].set_title(k)
    plt.savefig(outPath)
    # plt.show()

    print("Finished")

end = time.time()
print("Finished All")
print("Total time {}".format(end - start))

