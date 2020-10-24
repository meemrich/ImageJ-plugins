"""Rotates and crops a set of images from a directory and
saves the rotated cropped images in the input directory.
Designed for uniformly processing multiple exposures of a blots/gel
from a single imaging machine."""
#@ File[] (label="Select westerns", style="file") myImages

from ij import IJ, WindowManager
from ij.gui import GenericDialog, Roi, WaitForUserDialog
from ij.io import FileSaver
from ij.plugin.filter import Rotator
from ij.plugin.frame import RoiManager

from java.util import Date
from os import path
import re

class rotateCrop:
	def __init__(self):
		self.angle = None
		self.roi = int()

	def refProcessing(self, imp):
		"""process the refrence image and generate the roatation angle and roi used for cropping the subsequent images"""
		IJ.run("Rotate... ")
		self.angle = Rotator.getAngle()
		imp2 = imp.duplicate()
		wait = WaitForUserDialog("", "Please outline the membrane.\nPress 't' to add to the ROI manager. \nPress esc to cancel.")
		IJ.setTool("rectangle")
		rm = RoiManager.getInstance()
		if not rm:
			rm = RoiManager()
		rm.reset()
		wait.show()
		self.roi = rm.getRoi(0)
		imp.setRoi(self.roi)
		imp2 = imp.crop()
		return imp2

	def processImages(self, imp):
		"""rotate and crop images based on the angles and roi defined by the reference image"""
		IJ.run(imp, "Rotate... ", "angle="+str(self.angle)+ " grid=1 interpolation=Bilinear")
		imp.setRoi(self.roi)
		imp2 = imp.crop()
		return imp2

def saver(imp, output, name):
	"""save image plus object as tiff under the path output + name"""
	if path.exists(output) and path.isdir(output):
		filepath = path.join(output, name + ".tif")
		existTest = path.exists(filepath)
	else:
		print "Exiting.", output, "is not a valid directory"
	if existTest:
		gd = GenericDialog("results saver")
		text = "a file is already located at " + filepath + ". Please choose a new name"
		gd.addStringField(text, name + "-1")
		gd.showDialog()
		if gd.wasCanceled():
			print  "a file is already located at " + filepath + ". Renaming was cancelled by user."
			return
		newname = gd.getNextString()
		filepath = output + newname + ".tif"
	FileSaver(imp).saveAsTiff(filepath)
	print name, "saved successfully at ", str(filepath)

#open images and show all
myImagePaths = [str(i) for i in myImages]
imps = [IJ.openImage(i) for i in myImagePaths]
for imp in imps:
	imp.show()

output = IJ.getDirectory("image")
print "output directory:", output

#select a reference image. e.g. where the membrane is visible if a western,
picRef = WaitForUserDialog("", "Please select a reference blot.")
picRef.show()
imp = IJ.getImage()
imageName = re.split("\.\w{3}$", imp.getTitle())[0]
print "name of reference image:", imageName

rc = rotateCrop()
imp = rc.refProcessing(imp)
saver(imp, output, imageName + "_rotate_crop")

for imp in imps:
	name = re.split("\.\w{3}$", imp.getTitle())[0]
	if name != imageName:
		imp2 = rc.processImages(imp)
		saver(imp2, output, name + "_rotate_crop")
for imp in imps:
	imp.changes = False
	imp.close()

print "Image processing completed on", Date(), "."


