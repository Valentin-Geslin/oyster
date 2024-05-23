setImageType('OTHER');
runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2', '{"threshold":215,"requestedPixelSizeMicrons":20.0,"minAreaMicrons":1.0E6,"maxHoleAreaMicrons":1.0E9,"darkBackground":false,"smoothImage":true,"medianCleanup":true,"dilateBoundaries":true,"smoothCoordinates":true,"excludeOnBoundary":false,"singleAnnotation":true}')
resultingClass = getPathClass("TISSUE")
toChange = getAnnotationObjects().findAll{it.getPathClass() == null}
toChange.each{ it.setPathClass(resultingClass)}

import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathAnnotationObject
def project = getProject()
def imagesToExport = project.getImageList()
def separator = ","
def columnsToInclude = new String[]{""}
def exportType = PathAnnotationObject.class
def outputPath = System.getProperty("user.home")+"/OYSTER/results/tissue_area.csv"
def outputFile = new File(outputPath)
def exporter  = new MeasurementExporter()
                  .imageList(imagesToExport)            
                  .separator(separator)
                  .exportType(exportType)
                  .exportMeasurements(outputFile)    

print "Done!"