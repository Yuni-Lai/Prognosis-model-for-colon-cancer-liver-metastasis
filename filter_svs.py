# -*- coding: utf-8 -*-
import os
import fnmatch
import shutil
from collections import Counter

tilefolder='/data/backup/Yuni/CRC_Liver/05_tiles_output/'
slidefolder1='/data/backup/Yuni/CRC_Liver/08_new_SVS/'
slidefolder2='/data/backup/Yuni/CRC_Liver/liver_SVS/'
slidelist = os.listdir(slidefolder1)



newname=[slidename.replace(' ', '') for slidename in slidelist]
C1=Counter(newname)#578

fail=[]
for slidename in slidelist:
    try:
        os.rename(tilefolder+slidename,tilefolder+slidename.replace(' ', ''))
    except:
        print(f'{slidename} failed!')
        fail.append(slidename)
        continue

#remove file----------
for f in fail:
    print(tilefolder+f)
    shutil.rmtree(tilefolder+f)

for slidename in slidelist:
    if fnmatch.fnmatchcase(slidename, '*(1)*'):
        print(slidename)
        shutil.rmtree(tilefolder+slidename)
#        

for slidename in slidelist:
    if fnmatch.fnmatchcase(slidename, '*(4)*'):
        print(slidename)
        os.remove(tilefolder+slidename)






import os
import fnmatch
import shutil
from collections import Counter

slidefolder='/data/backup2/Yuni/CRC_Liver/04_multi_regions/CRC_LM_10/'
slidelist = os.listdir(slidefolder)

def get_image_names(folder):
    imagefolder=os.listdir(folder)
    image_names=[]  
    for f in imagefolder:
        images=os.listdir(folder+f)
        images = [f+'/'+im for im in images]
        image_names.extend(images)
    return image_names

slidelist=get_image_names(slidefolder)

fail=[]
for slidename in slidelist:
    try:
        os.rename(slidefolder+slidename,slidefolder+slidename.replace(' ', '_'))
    except:
        print(f'{slidename} failed!')
        fail.append(slidename)
        continue     



#-----remove empty files-------------------------------------------
#'/data/backup/Yuni/CRC_Liver/01_tiles_output/04_tiles_output/'
#'/data/backup/Yuni/CRC_Liver/01_tiles_output/05_tiles_output/' 
#/data/backup2/Yuni/CRC_Liver/01_tiles_output/04_tiles_output/
#/data/backup2/Yuni/CRC_Liver/01_tiles_output/05_tiles_output/
#have a test first!!!!
tilefolder='/data/backup2/Yuni/CRC_Liver/01_tiles_output/05_tiles_output/'
#tilelist = get_image_names(tilefolder)
tilelist = os.listdir(tilefolder)
for tile in tilelist:
    if fnmatch.fnmatchcase(tile, '*.pkl') and os.path.getsize(tilefolder+tile)>0:
        pass
    else:
        print(f'remove: {tile}')

#remove it if sure!!!

for tile in tilelist:
    #tile=1026806.pkl
    if fnmatch.fnmatchcase(tile, '*.pkl'):
        if os.path.getsize(tilefolder+tile)==0:
            print(f'remove: {tile}')
            os.remove(tilefolder+tile)
            
            
#    else:
#        print(f'remove: {tile}')
#        shutil.rmtree(tilefolder+tile)

