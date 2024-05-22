# Detection, quantification and characterisation of bacterial infection by _Vibrio aestuarianus_, stained by immunohistochemistry, in the pacific oyster _Magallana gigas_ by digital images analysis method

This repository contains the source code for the paper "Detection, quantification and characterisation of bacterial infection by _Vibrio aestuarianus_, stained by immunohistochemistry, in the pacific oyster _Magallana gigas_ by digital images analysis method".
Published version: [https://.../](https://.../)

Please cite the [article](https://.../) if you use this code or build on top of it.

## Requirements
- QuPath (0.4)
- Python (3.8)
- Fiji (2.15)
- R (4.3.1)
- HistoQC (https://github.com/choosehappy/HistoQC)

## How to use the code
- Clone the repository locally (do not modify file/folder structure within the repository)
- Download digital slides to the WSI folder (https://seanoe.org/data/00501/61299/)
- Open QuPath project and load the WSI to the project
- Run "export_metadata.groovy" for the entire project with QuPath and copy the results to /results/image_metadata.txt
- Run "tissue_detection.groovy" for the entire project with QuPath
- Run "export_tile.groovy" for the entire project with QuPath
- Run "ihc_detection_tile.py" with Python 
- Run HistoQC on WSI files (use "light" configuration)
- Run "sort_histoqc_output.py" with Python
- Run "crop_image.py" with Python
- Run "ihc_detection_image.py" with Python
- Run "colocalization.ijm" with Fiji
- Run "dispersion.ijm" with Fiji
- Run "data_analysis.R" with R
