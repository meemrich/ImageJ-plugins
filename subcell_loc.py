#1. reads in czi files with the nucleus stained in one channel and whole cells visible in another
#2. measures intensity and area of cells, nuclei, and cytosol in the channels of interest
#3. outputs measurements as a csv file

#@ File (label="Select the input file", description="input .czi file location") imagefile

from ij import IJ, ImagePlus, Prefs, WindowManager
from ij.process import ImageStatistics as IS
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import Overlay, Roi, ShapeRoi, WaitForUserDialog, GenericDialog
# read in and display ImagePlus object(s)
from loci.plugins import BF
from loci.common import Region
from loci.plugins.in import ImporterOptions
from loci.plugins.util import LociPrefs
from loci.formats import ImageReader, MetadataTools
# ZEISS STUFF
from loci.formats.in import ZeissCZIReader, DynamicMetadataOptions
#regular expressions
import re
#saving
from os import path

import sys
import collections

#file management
from java.io import File
from java.lang import System
from java.text import SimpleDateFormat

#detect if run via Gui
from java.awt import GraphicsEnvironment


def binarize(imp, channel):
	"""convert image plus to binary imp"""
	ip = imp.getProcessor()
	stats = IS.getStatistics(ip, IS.MIN_MAX, imp.getCalibration())
	maxVal = stats.max
	imp.setDisplayRange(450, maxVal)
	IJ.run(imp, "Apply LUT", "")
	IJ.run(imp, "Gaussian Blur...", "sigma=1")
	IJ.run(imp, "Make Binary", "BlackBackground")
	return imp

def channelSelector(channelCount, channelIDs = None):
	"""identify which channels contain which signals in a given image and name the labels for each"""
	if not channelIDs:
		gd = GenericDialog("channel selector")
		gd.addMessage("Please input the channel names and indicies identifying the nucleus, the cell, and the protein(s) of interest.\nIf there is more than one protein of interest, please add each as comma separated values.\nInput data in the format 'label = number' as in the defaults.")
		gd.addStringField("nuclear channel", 'DAPI = 1')
		gd.addStringField("cell channel", 'Cherry = 3')
		gd.addStringField("protein(s) of interest", 'GFP = 2, Cherry = 3')
		gd.showDialog()
		if gd.wasCanceled():
			IJ.log("Channel selection was cancelled. Exiting")
			IJ.error("Channel selection was cancelled. Exiting")
			return
		channelIDs = [gd.getNextString() for i in range(3)]
		proteinsOfinterest = str(channelIDs.pop()).split(",")
		channelIDs += proteinsOfinterest
		channelDict = {}
	idList = ['nucleus', 'cell'] + [i for i in range(len(channelIDs) - 2)]
	while channelIDs:
		val = idList.pop()
		protein = str(channelIDs.pop()).split("=")
		try:
			channelNo = int(protein[1])
		except:
			IJ.log("Non-numeric channel input. Exiting")
			IJ.error("Non-numeric channel input. Exiting")
			return
		if 1 <= channelNo <= channelCount:
			channelDict[val] = [protein[0], channelNo]
		else:
			IJ.log("Channel input outside bounds. Exiting")
			IJ.error("Channel input outside bounds. Exiting")
			return
	#print channelDict
	return channelDict

def getCZIinfo(imagefile):
	"""import czi info incl. image dimensions and series length"""
	CZIinfo = {}
	options = DynamicMetadataOptions()
	options.setBoolean("zeissczi.autostitch", False)
	options.setBoolean("zeissczi.attachments", False)
	czireader = ZeissCZIReader()
	czireader.setFlattenedResolutions(False)
	czireader.setMetadataOptions(options)
	czireader.setId(imagefile)
	CZIinfo['seriesCount'] = czireader.getSeriesCount()
	CZIinfo['SizeC'] = czireader.getSizeC()
	CZIinfo['SizeX'] = czireader.getSizeX()
	CZIinfo['SizeY'] = czireader.getSizeY()
	czireader.close()
	#print "CZIinfo generated"
	return CZIinfo

def imageOpener(imagefile):
	# Set the preferences in the ImageJ plugin
	Prefs.set("bioformats.zeissczi.allow.autostitch", str(False).lower())
	Prefs.set("bioformats.zeissczi.include.attachments", str(True).lower())
	# read in and display ImagePlus(es) with arguments
	options = ImporterOptions()
	options.setOpenAllSeries(True)
	options.setConcatenate(False)
	options.setId(imagefile)
	imps = BF.openImagePlus(options)
	return imps

def processimagefile(imagefile):
	"""function to take in input file and process it. Runs the other functions in this plugin."""
	sep = System.getProperty("file.separator")
	output = imagefile.getParent() + sep
	imagefile = str(imagefile)
	splitFile = imagefile.split("/")
	name = splitFile[-1].split(".")[0]
	#process image metadata and open image file
	CZIinfo = getCZIinfo(imagefile)
	imp = imageOpener(imagefile)
	outputArray = []
	#If image has fewer than 2 channels exit
	#otherwise the user identifies which channels label cellular compartments and which label proteins of interest
	if CZIinfo['SizeC'] < 2:
		IJ.log("A minimum of 2 channels is required for identification of cells and nuclei. Exiting")
		IJ.error("A minimum of 2 channels is required for identification of cells and nuclei. Exiting")
		return
	else:
		channels = channelSelector(CZIinfo['SizeC'])
	if not channels:
		return
	#loop through all frames in the image. In each frame, identify all cells and
	#their respective nuclei and output an array containing the area and intensity of each
	for frame in range(CZIinfo['seriesCount']):
		outputArray.append(imageProperties().imageProcessor(imp[frame], CZIinfo['SizeX'], CZIinfo['SizeY'], channels))
		IJ.log("finished processing {0} out of {1} frames.".format(frame + 1, CZIinfo['seriesCount']))
	IJ.log("image processing finished")
	log = IJ.getLog()
	resultsSaver(log, output, name + "_log", ".txt")
	table = resultsTablemaker(outputArray)
	if ResultsTable.size(table) == 0:
		IJ.log("No transfected cells found. Bye.")
		return
	resultsSaver(table, output, name + "_results", ".csv")
	WindowManager.closeAllWindows()

def resultsTablemaker(outputArray):
	"""convert flattened list to results table"""
	IJ.log("filling results table")
	table = ResultsTable()
	colnames = collections.deque(["Name", "snapNo", "Channel", "ROI", "Area", "Mean", "Median"])
	colLen = len(colnames)
	for image in outputArray:
		if image is None:
			continue
		for cell in image:
			for i, value in enumerate(cell):
				if i%colLen == 0:
					table.incrementCounter()
				col = colnames.popleft()
				table.addValue(col, value)
				colnames.append(col) 
	return table

def resultsSaver(item, output, name, extension):
	"""save results table and log files"""
	filepath = output + name + extension
	level = 0
	while File(filepath).exists():
		level += 1
		filepath = output + name + str(level) + extension
	if extension == ".csv":
		item.save(filepath)
	elif extension == ".txt":
		with open(filepath, 'w') as logSaver:
			logSaver.write(item)
	IJ.log("Output saved to {0}".format(filepath))

class imageProperties:
	"""set of functions to identify the cells and nuclei in a given imp and measure the intensity of these relative to bg"""
	def __init__(self):
		self.images = {}
		self.imageLabels = {}
		self.cellOverlay = Overlay()
		self.DAPIoverlay = Overlay()
	
	def cytoplasmMaker(self, cellCount):
		finalOverlay = Overlay()
		for i in range(cellCount):
			roi = self.cellOverlay.get(i)
			Rname = Roi.getName(roi)
			tempNo = re.split("_", Rname)[-1]
			roi2 = self.DAPIoverlay.get("nucleus_" + tempNo)
			cytoplasm = ShapeRoi(roi).not(ShapeRoi(roi2))
			Roi.setName(cytoplasm, "cytoplasm_" + tempNo)
			finalOverlay.add(roi)
			finalOverlay.add(roi2)
			finalOverlay.add(cytoplasm)
		if Overlay.size(finalOverlay) == 0:
			IJ.log("no transfected cells identified in image {0} frame {1} after overlap filtering".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None
		finalCells = Overlay.size(finalOverlay)/3
		if isinstance(finalCells, (int, long)):
			IJ.log("{0} cells identified in image {1} frame {2} after filtering".format(finalCells, self.imageLabels['snapName'], self.imageLabels['snapNo']))
		else:
			IJ.log("incorrect cell:nucleus:cytoplasm ratio in image {0} frame {1} after mask filtering".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None
		return finalOverlay
	
	def imagePlusmaker(self, imp, channels):
		"""fill the self.imageLabels and self.images dictionaries with image data and the imagePlus' of the cell and nucleus channels.
		Need to run this to create the inherited variables used in the other functions"""
		self.imageLabels['snapName'] = imp.getTitle().split(".")[0]
		if "#" in imp.getTitle():
			self.imageLabels['snapNo'] = imp.getTitle().split("#")[1]
		else:
			IJ.log("Unable to identify snap number, using entire image name instead.")
			self.imageLabels['snapNo'] = imp.getTitle()
		IJ.log("beginning processing of image {0} frame {1}".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
		nucleusLabel = channels['nucleus'][0]
		cellLabel = channels['cell'][0]
		self.images[nucleusLabel] = ImagePlus(nucleusLabel, imp.getImageStack().getProcessor(channels['nucleus'][1])).duplicate()
		self.images[cellLabel] = ImagePlus(cellLabel, imp.getImageStack().getProcessor(channels['cell'][1])).duplicate()
		return imp
	
	def itemID(self, imp, item, shape):
		"""manually identify objects and add them to an overlay"""
		if imp.getOverlay():
			imp.getOverlay().clear()
		rm = RoiManager.getInstance()
		if not rm:
			rm = RoiManager()
		rm.reset()
		imp.show()
		wait = WaitForUserDialog("", "Please draw "+ item + ".\nPress 't' to add to the ROI manager.")
		if shape == "r":
			IJ.setTool("rectangle")
		else:
			IJ.setTool("polygon")
		wait.show()
		imp.hide()
		if wait.escPressed():
			IJ.log("{0} drawing cancelled for image {1} frame {2}".format(item, self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return
		elif rm.getCount() == 0:
			IJ.log("no {0} drawn for image {1} frame {2}".format(item, self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return
		else:
			IJ.log("{0} drawn for image {1} frame {2}".format(item, self.imageLabels['snapName'], self.imageLabels['snapNo']))
			rm.moveRoisToOverlay(imp)
			rm.reset()
			imp.hide()
	
	def measureImage(self, imp, overlay, channels):
		overlaySize = Overlay.size(overlay)
		imp.setOverlay(overlay)
		table = []
		for label, channel in channels.items():
			if label == 'nucleus' or label == 'cell':
				continue
			imp.setSlice(channel[1])
			temp = []
			for i in range(overlaySize):
				roi = overlay.get(i)
				roi.setImage(imp)
				Rname = Roi.getName(roi)
				roiStat = roi.getStatistics()
				temp.append(self.imageLabels['snapName'])
				temp.append(self.imageLabels['snapNo'])
				temp.append(channel[0])
				temp.append(Rname)
				temp.append(roiStat.area)
				temp.append(roiStat.mean)
				temp.append(roiStat.median)
			table.append(temp)
		return table
	
	def nucleiFilter(self, cellPoints, cellCount, nucleiCount):
		"""iterate across all nuclei and delete those that don't overlap ChREBP-transfected cells"""
		if cellCount < nucleiCount:
			i = 0
			while self.DAPIoverlay.contains(self.DAPIoverlay.get(i)):
				roi = self.DAPIoverlay.get(i)
				nucPoints = set(roi.getContainedPoints())
				if len(cellPoints & nucPoints) < 0.9 * len(nucPoints):
					self.DAPIoverlay.remove(roi)
					nucleiCount -= 1
					if nucleiCount == 0:
						IJ.log("no nuclei found to be coincident with called transfected cells in image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
						self.images[cellLabel].show()
						self.itemID(self.images[nucleusLabel], "nuclei", "p")
						self.images[cellLabel].hide()
						self.DAPIoverlay = self.images[nucleusLabel].getOverlay()
						nucleiCount = Overlay.size(self.DAPIoverlay)
				else:
					i += 1
		return nucleiCount
	
	def overlayMaker(self, cellLabel, nucleusLabel, width, height):
		self.itemID(self.images[cellLabel], "cells", "p")
		self.cellOverlay = self.images[cellLabel].getOverlay()
		if not self.cellOverlay:
			IJ.log("no cells found in image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None , None, None
		cellCount = Overlay.size(self.cellOverlay)
		cellRoi = ShapeRoi(self.cellOverlay.get(0))
		for i in range(1, cellCount):
			shape = ShapeRoi(self.cellOverlay.get(i))
			cellRoi = cellRoi.or(shape)
		cellPoints = set(cellRoi)
		notCell = ShapeRoi(Roi(0, 0, width, height)).xor(cellRoi)
		ip = self.images[nucleusLabel].getProcessor()
		ip.setValue(0)
		ip.setRoi(notCell)
		ip.fill(ip.getMask())
		binarize(self.images[nucleusLabel], "DAPI")
		"""create an overlay of _all_ nuclei"""
		IJ.run(self.images[nucleusLabel], "Analyze Particles...", "size=1500-10000 circularity=0.4-1.00 show=Overlay include")
		self.DAPIoverlay = self.images[nucleusLabel].getOverlay()
		if not self.DAPIoverlay:
			IJ.log("no nuclei automatically called found in image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			self.images[cellLabel].show()
			self.itemID(self.images[nucleusLabel], "nuclei", "p")
			self.images[cellLabel].hide()
			self.DAPIoverlay = self.images[nucleusLabel].getOverlay()
		if not self.DAPIoverlay:
			IJ.log("no nuclei found in image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None , None, None
		nucleiCount = Overlay.size(self.DAPIoverlay)
		IJ.log("initial cells: {0}; initial nuclei: {1}".format(cellCount, nucleiCount))
		cellRoi = ShapeRoi(self.cellOverlay.get(0))
		return cellPoints, cellCount, nucleiCount
		
	def roiFilter(self, cellCount, nucleiCount, outerType):
		"""to identify which nucleus belongs to which cell or vice verse in cases of more than one possible nucleus per cell"""
		if outerType == "cell":
			outerOverlay, innerOverlay, outerCount, innerCount = self.cellOverlay, self.DAPIoverlay, cellCount, nucleiCount
		else:
			outerOverlay, innerOverlay, outerCount, innerCount = self.DAPIoverlay, self.cellOverlay, nucleiCount, cellCount
		i = 0
		while i < outerCount:
			overlapIndex, overlapTestList = [], []
			roi = outerOverlay.get(i)
			outerPoints = roi.getContainedPoints()
			j = 0
			while j < innerCount:
				roi2 = innerOverlay.get(j)
				innerPoints = roi2.getContainedPoints()
				overlapTest = len(set(outerPoints) & set(innerPoints))
				if overlapTest != 0:
					overlapIndex.append(j)
					overlapTestList.append(overlapTest)
				if j == innerCount - 1 and len(overlapIndex) == 1:
					"""if a given cell/nucleus ovelaps exactly one nucleus/cell- label the respective cell & nucleus with matching
					identifiers starting with the str unique"""
					roi4 = innerOverlay.get(overlapIndex[0])
					roi4.setName("unique" + outerType + str(i)+ str(j))
					roi.setName("unique" + outerType + str(i)+ str(j))
				if j == innerCount - 1 and len(overlapIndex) > 1:
					"""if a given cell/nucleus ovelaps more than one nucleus/cell- label the respective cell-nucleus pair
					with the greatest overlap with matching identifiers"""
					maxInd = overlapTestList.index(max(overlapTestList))
					tempMax = overlapIndex[maxInd]
					roi4 = innerOverlay.get(tempMax)
					if roi4.getName() is None:
						roi4.setName(outerType + str(i)+ str(j))
						roi.setName(outerType + str(i)+ str(j))
					else:
						"""if the labelled cell-nucleus pair happen to contain a unique cell-nucleus pair, it is unlikely to be a true
						pair so assign matching identifiers to the cell-nucleus pair with the next smallest overlapping area"""
						while "unique" in roi4.getName():
							overlapIndex.remove(tempMax)
							overlapTestList.remove(max(overlapTestList))
							if len(overlapIndex) == 0:
								break
							maxInd = overlapTestList.index(max(overlapTestList))
							tempMax = overlapIndex[maxInd]
							roi4 = innerOverlay.get(tempMax)
						if len(overlapIndex) > 0:
							roi4.setName(outerType + str(i)+ str(j))
							roi.setName(outerType + str(i)+ str(j))
				j += 1
			if len(overlapIndex) == 0:
				outerOverlay.remove(roi)
				outerCount -= 1
			else:
				if i == outerCount - 1:
					for k in range(innerCount):
						roi3 = innerOverlay.get(i)
						if innerOverlay.contains(roi3) == True:
							if roi3.getName() is None:
								innerOverlay.remove(roi3)
								innerCount -= 1
								if innerOverlay == 0:
									IJ.log("no single transfected cells overlapped a single transfected nucleus in image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
									return None, None
						else:
							innerCount = k
							break
				i += 1
		if outerType == "cell":
			self.cellOverlay, self.DAPIoverlay, cellCount, nucleiCount = outerOverlay, innerOverlay, outerCount, innerCount
		else:
			self.DAPIoverlay, self.cellOverlay, nucleiCount, cellCount = outerOverlay, innerOverlay, outerCount, innerCount
		return cellCount, nucleiCount
	
	def unmatchedFilter(self, cellCount, nucleiCount):
		"""to remove all cells that do not overlap with an precisely one nucleus and vice versa"""
		for i in range(cellCount):
			roi = self.cellOverlay.get(i)
			Rname = roi.getName()
			if self.DAPIoverlay.contains(self.DAPIoverlay.get(Rname)) == False:
				roi.setName("unmatched" + str(i))
		for i in range(nucleiCount):
			roi = self.DAPIoverlay.get(i)
			Rname = roi.getName()
			if self.cellOverlay.contains(self.cellOverlay.get(Rname)) == False:
				roi.setName("unmatched" + str(i))
		for i in range(cellCount):
			roi = self.cellOverlay.get(i)
			Rname = roi.getName()
			if "unmatched" in Rname:
				for j in range(nucleiCount):
					roi2 = self.DAPIoverlay.get(i)
					if self.DAPIoverlay.contains(roi2) == True:
						Rname2 = roi2.getName()
						if "unmatched" in Rname2:
							cellPoints = roi.getContainedPoints()
							nucPoints = roi2.getContainedPoints()
							overlapTest = bool(set(cellPoints) & set(nucPoints))
							if overlapTest == True:
								roi.setName("name" + str(i))
								roi2.setName("name" + str(i))
		i = 0
		while i < cellCount:
			roi = self.cellOverlay.get(i)
			Rname = roi.getName()
			roi2 = self.DAPIoverlay.get(Rname)
			if self.DAPIoverlay.contains(roi2) == True:
				roi.setName("cell_" + str(i))
				roi2.setName("nucleus_" + str(i))
				i += 1
			else:
				self.cellOverlay.remove(roi)
				cellCount -= 1
		return cellCount
	
	def imageProcessor(self, imp, width, height, channels):
		"""function to run the other functions in this class and output an array of image measurements"""
		self.imagePlusmaker(imp, channels)
		nucleusLabel = channels['nucleus'][0]
		cellLabel = channels['cell'][0]
		cellPoints, cellCount, nucleiCount = self.overlayMaker(cellLabel, nucleusLabel, width, height)
		if cellPoints is None:
			return None
		nucleiCount = self.nucleiFilter(cellPoints, cellCount, nucleiCount)
		if nucleiCount == 0:
			IJ.log("no nuclei found in image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None
		cellCount, nucleiCount = self.roiFilter(cellCount, nucleiCount, "cell")
		if cellCount == None:
			return None
		cellCount, nucleiCount = self.roiFilter(cellCount, nucleiCount, "DAPI")
		if cellCount == None:
			return None
		cellCount = self.unmatchedFilter(cellCount, nucleiCount)
		finalOverlay = self.cytoplasmMaker(cellCount)
		if finalOverlay == None:
			return None
		self.itemID(self.images[cellLabel], "bg", "r")
		bg = self.images[cellLabel].getOverlay().get(0)
		if not bg:
			IJ.log("no bg drawn for image {0} frame {1}.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None
		bg.setName("bg")
		Overlay.add(finalOverlay, bg)
		table = self.measureImage(imp, finalOverlay, channels)
		if len(table) > 0:
			return table
		else:
			IJ.log("table not generated for image {0} frame {1} even though transfected cells overlapped specific nuclei.".format(self.imageLabels['snapName'], self.imageLabels['snapNo']))
			return None

if __name__ in ["__builtin__", "__main__"]:
	processimagefile(imagefile)
