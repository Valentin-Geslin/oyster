from glob import glob
import os
import shutil

curr_dir = os.getcwd()
thumbs_dir = curr_dir+os.path.sep+"images_histoqc"+os.path.sep
first_masks_dir = curr_dir + os.path.sep + "mask_histoqc" + os.path.sep
masks_dir = first_masks_dir + os.path.sep + "mask_use" + os.path.sep

for fold in [thumbs_dir,first_masks_dir,masks_dir]:
    if not os.path.isdir(fold): os.mkdir(fold)

histoqc_dir = curr_dir.replace("OYSTER", "HistoQC")

subfolders = [(f.path, os.path.getmtime(f.path)) for f in os.scandir(histoqc_dir) if f.is_dir() and "histoqc_output_" in f.path]
assert len(subfolders)>0, "There are no directories with HistoQC results"

results_dir = subfolders
if len(subfolders)>1: results_dir = list(map(max,*zip(subfolders)))
for subfold in glob(results_dir[0][0]+os.path.sep+"*"):
    if ".ndpi" in subfold:
        os.rename(subfold, subfold.replace(".ndpi", ""))

    for image in glob(subfold+os.path.sep+"*.png"):
        if ".ndpi_thumb" in image and "small" not in image:
            shutil.copyfile(image, thumbs_dir+image.split(os.path.sep)[-1].replace(".ndpi_thumb", "")
                            .replace(".png","").strip()+".png")
        if ".ndpi_mask_use" in image:
            shutil.copyfile(image, masks_dir + image.split(os.path.sep)[-1].replace(".ndpi_mask_use", "")
                            .replace(".png","").strip()+".png")