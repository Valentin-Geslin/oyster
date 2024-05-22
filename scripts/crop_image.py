# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 13:29:16 2022

@author: ValentinGeslin
"""

import cv2
import glob
import numpy as np
import os

curr_dir = os.getcwd()
img_path = curr_dir+os.path.sep+"images_histoqc"+os.path.sep
img_masked_path = curr_dir + os.path.sep + "mask_histoqc" + os.path.sep
for file in glob.glob(img_path+'/*.png'):
    img = cv2.imread(file)
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    lower = np.array([132, 30, 0])
    upper = np.array([179, 255, 255])
    mask = cv2.inRange(hsv, lower, upper)
    _, mask = cv2.threshold(mask, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
    img_masked = cv2.bitwise_and(img, img, mask=mask)
    basename = os.path.basename(file)
    img_name = os.path.splitext(basename)[0]
    cv2.imwrite(os.path.join(img_masked_path, img_name+".png"), img_masked)