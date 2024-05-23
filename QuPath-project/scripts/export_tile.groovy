def imageData = getCurrentImageData()
def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def pathOutput = buildFilePath(System.getProperty("user.home")+'/OYSTER/tiles/', name)
mkdirs(pathOutput)
double requestedPixelSize = 1
double pixelSize = imageData.getServer().getPixelCalibration().getAveragedPixelSize()
double downsample = requestedPixelSize / pixelSize
new TileExporter(imageData)
    .downsample(downsample)   
    .imageExtension('.png')   
    .tileSize(400) 
    .annotatedTilesOnly(true)
    .overlap(0)
    .writeTiles(pathOutput)
print 'Done!'
