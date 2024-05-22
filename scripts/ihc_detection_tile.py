# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:00:41 2023

@author: Valentin_Geslin
"""

import csv
import cv2
import glob
import numpy as np
import os

curr_dir = os.getcwd()
tile_path = curr_dir+"/tiles/"
tile_masked_path = curr_dir+"/tiles_masked/"

for image in os.listdir(tile_path):
    path = os.path.join(tile_path, image)
    masked_path = os.path.join(tile_masked_path, image)
    if not os.path.isdir(masked_path): os.makedirs(masked_path)
    with open(curr_dir+"/results/IHC_tiles/mask_area_{}.csv".format(image), "a", newline="") as csvfile:
        fieldnames = ["image", "tile_name", "total_pixels", "black_pixels", "non_black_pixels"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for file in glob.glob(path + '/*.png'):
            tile_img = cv2.imread(file)
            hsv = cv2.cvtColor(tile_img, cv2.COLOR_BGR2HSV)
            lower = np.array([132, 30, 0])
            upper = np.array([179, 255, 255])
            mask = cv2.inRange(hsv, lower, upper)
            _, mask = cv2.threshold(mask, thresh=180, maxval=255, type=cv2.THRESH_BINARY)
            tile_masked = cv2.bitwise_and(tile_img, tile_img, mask=mask)
            basename = os.path.basename(file)
            tile_name = os.path.splitext(basename)[0]
            cv2.imwrite(os.path.join(masked_path, tile_name + "_masked.png"), tile_masked)
            height, length, dim = tile_masked.shape
            nb_pix = height * length * dim  # total number of pixels
            nb_black_pix = np.count_nonzero(tile_masked == 0)
            nb_non_black_pix = np.count_nonzero(tile_masked != 0)
            writer.writerow({"image": image, "tile_name": tile_name, "total_pixels": nb_pix, "black_pixels": nb_black_pix,
                             "non_black_pixels": nb_non_black_pix})