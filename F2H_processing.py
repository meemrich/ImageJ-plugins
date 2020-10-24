#### to quantify F2H czi images
### automatically identifies nucleus and array and outputs the area, mean and median intensity of each feature identified to a csv file
### the ROIs for each respective image are saved to the output directory to allow confirmation of accurate feature calling

#@ File[] (label="Select the input directories", style="directories") inputDir

#ImageJ stuff
from ij import IJ, ImagePlus, Prefs, WindowManager
from ij.process import ImageStatistics as IS
from ij.io import FileSaver
from ij.plugin import ImageCalculator, filter
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import Overlay, Roi, ShapeRoi, GenericDialog

# czi and image import stuff
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats.in import ZeissCZIReader, DynamicMetadataOptions

#python stuff (regular expressions etc.)
import re
import collections

#file management
from java.io import File
from java.lang import System
from java.text import SimpleDateFormat

#GUI stuff
from java.awt import GridBagLayout, GridBagConstraints
from javax.swing import JDialog, JFrame, JPanel, JLabel, JTextField, BorderFactory, JButton

class frameMaker():
	def __init__(self):
		self.Imageinfo = {}
		self.dialog = JDialog()
		self.panel = JPanel(border = BorderFactory.createEmptyBorder(10,10,10,10))
		self.out = {}
		
	def okayPressed(self, event):
		length = len(self.panel.getComponents())
		upperFinder = re.compile("^upper.*")
		for i in range(1, length):
			item = self.panel.getComponents()[i]
			if not isinstance(item, JTextField):
				continue
			label = self.panel.getComponents()[i-1].getText() #labels are located one index behind their respective textboxes in the panel
			if label in ["nuclear", "prey", "bait"]:
				bounds = [1, self.Imageinfo['SizeC']]
			elif label in ["lower nuclear area", "upper nuclear area"]:
				bounds = [1, self.Imageinfo['SizeX']*self.Imageinfo['SizeY']]
			elif label in ["lower array area", "upper array area"]:
				bounds = [1, self.out['lower nuclear area']]
			elif label in ["lower thresh", "upper thresh"]:
				bounds = [0, 255]
			try:
				val = int(item.getText())
			except:
				IJ.log("Non-numeric input for {0}. Please input an integer between {1} and {2}.".format(label, bounds[0], bounds[1]))
				IJ.error("Non-numeric input for {0}. Please input an integer between {1} and {2}.".format(label, bounds[0], bounds[1]))
				return
			if upperFinder.match(label): #check if label begins with upper
				lowerVal = int(self.panel.getComponents()[i-2].getText())
				if val <= lowerVal:#the upper value of any pair of upper and lower values cannot be less than or equal to the lower value
					sectionLabel = label.split("upper ")[1]
					IJ.log("Input value outside permitted range. Upper {0} must be greater than lower {0}.".format(sectionLabel))
					IJ.error("Input value outside permitted range. Upper {0} must be greater than lower {0}.".format(sectionLabel))
					return
			if bounds[0] <= val <= bounds[1]:
				self.out[label] = val
			else:
				IJ.log("Input value {0} outside of range. Please input an integer between {1} and {2} in {3}".format(val, bounds[0], bounds[1], label))
				IJ.error("Input value {0} outside of range. Please input an integer between {1} and {2} in {3}".format(val,bounds[0], bounds[1], label))
				return
		print self.out, "okaypressed"
		self.dialog.dispose()

	def cancelPressed(self, event):
		IJ.error("Parameter selection cancelled. Exiting.")
		self.dialog.dispose()
	
	def buttonMaker(self, gb):
		gbc = GridBagConstraints()
		okayButton = JButton("okay", actionPerformed = self.okayPressed)
		gbc.gridx = 0
		gbc.anchor = GridBagConstraints.WEST
		gb.setConstraints(okayButton, gbc)
		self.panel.add(okayButton)
		cancelButton = JButton("cancel", actionPerformed = self.cancelPressed)
		gbc.gridx = 1
		gb.setConstraints(cancelButton, gbc)
		self.panel.add(cancelButton)
	
	def inputMaker(self, gb, vals):
		"""add the labels and text input fields for the respective sections"""
		gbc = GridBagConstraints()
		for val in vals:
			gbc.anchor = GridBagConstraints.WEST
			gbc.gridx = 0 
			label = JLabel(val[0])
			gb.setConstraints(label, gbc)
			self.panel.add(label)
			gbc.gridx = 1
			gb.setConstraints(val[1], gbc)
			self.panel.add(val[1])
	
	def sectionLabel(self, gb, text):
		"""create labels 2 columns wide to separate types of input"""
		gbc = GridBagConstraints()
		gbc.gridx = 0
		gbc.anchor = GridBagConstraints.WEST
		gbc.gridwidth = 2
		sectionLabel = JLabel(text)
		gb.setConstraints(sectionLabel, gbc)
		self.panel.add(sectionLabel)
	
	def idBuilder(self, IDs):
		IDs['channel'] = [("nuclear", JTextField("3", 5)), ("prey", JTextField("2", 5)), ("bait", JTextField("1", 5)), "<html>Input the indices of the indicated channels.</html>"]
		IDs['nuc'] = [("lower nuclear area", JTextField("4000", 5)), ("upper nuclear area", JTextField("20000", 5)), "<html> <br/>Input the minimum and maximum nuclear area <br/> (in pixels) for nucleus calling.</html>"]
		IDs['array'] = [("lower array area", JTextField("5", 5)), ("upper array area", JTextField("200", 5)), "<html> <br/>Input the minimum and maximum array area <br/> (in pixels) for array calling.</html>"]
		IDs['thresh'] = [("lower threshold", JTextField("97", 5)), ("upper threshold", JTextField("195", 5)), "<html> <br/>Input the values (between 0 and 255) to threshold <br/> the bait images for array calling.</html>"]
		return IDs
		
	def dialogBuilder(self, ImageInfo):
		self.Imageinfo = ImageInfo
		"""set up basic panel with slight border and establish the gb layout"""
		gb = GridBagLayout()
		self.panel.setLayout(gb)
		"""fill the self.IDs hash table of labels and input boxes"""
		IDs = self.idBuilder({})
		"""as dictionaries are not ordered in python 2, create separate list of keys to allow looping through self.IDs in the desired order"""
		idList = ['channel', 'nuc', 'array', 'thresh']
		"""use idList to loop through self.IDs and add an instruction label and labelled text boxes for each section"""
		for item in idList:
			temp = IDs[item]
			instructions = temp.pop()
			self.sectionLabel(gb, instructions)
			self.inputMaker(gb, temp)
		"""add the okay and cancel buttons for the dialog box"""
		self.buttonMaker(gb)
		"""setup the frame"""
		self.dialog.add(self.panel)
		self.dialog.pack()
		self.dialog.setLocationRelativeTo(None)
		self.dialog.setModal(True)
		self.dialog.setVisible(True)
		return self.out
		
def CZIopener(imagefile):
	"""import czi info incl. image dimensions and series length"""
	options = ImporterOptions()
	options.setAutoscale(True)
	options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
	options.setId(imagefile)
	imp = BF.openImagePlus(options)[0]
	return imp

def findnucleus(imp, parameters):
	"""call nuclei"""
	IJ.run(imp, "Gaussian Blur...", "sigma=3")
	IJ.run(imp, "Enhance Contrast", "saturated=0.35")
	IJ.run(imp, "Apply LUT", "")
	IJ.run(imp, "Auto Threshold", "method=Default white")
	IJ.run(imp, "Make Binary", "BlackBackground")
	IJ.run(imp, "Analyze Particles...", "size=" + str(parameters['lower nuclear area']) + "-" + str(parameters['upper nuclear area']) + " circularity=0.50-1.00 show=Overlay include")
	DAPIoverlay = imp.getOverlay()
	return DAPIoverlay

def findarray(images, DAPIoverlay, totalnuclei, parameters):
	"""use the nuclear channel to make a mask of all non-nuclear areas in the image, then identify the arrays"""
	nucShape = ShapeRoi(DAPIoverlay.get(0))
	i = 1
	while i < totalnuclei:
		shape = ShapeRoi(DAPIoverlay.get(i))
		nucShape = nucShape.or(shape)
		i += 1
	bgShape = ShapeRoi(Roi(0, 0, images['bait'].getWidth(), images['bait'].getHeight()))
	bgShape = bgShape.not(nucShape)
	ip = images['bait'].getProcessor()
	ip.setValue(0)
	ip.setRoi(bgShape)
	ip.fill(ip.getMask())
	IJ.run(images['bait'], "8-bit", "")
	IJ.setThreshold(images['bait'], parameters['lower threshold'], parameters['upper threshold'], "Black & White")
	IJ.run(images['bait'], "Analyze Particles...", "size=" + str(parameters['lower array area']) + "-" + str(parameters['upper array area']) + "pixel circularity=0.50-1.00 show=Overlay include")
	baitoverlay = images['bait'].getOverlay()
	return baitoverlay

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
	return CZIinfo

def imageProcessor(imagefile, parameters, outputDir):
	imp = CZIopener(imagefile)
	images, imageLabels = {}, {}
	imageLabels['snapName'] = imp.getTitle().split(".")[0]
	imageLabels['imagefile'] = imagefile
	imageLabels['date'] = SimpleDateFormat("yyyy/MM/dd").format(File(imagefile).lastModified())
	IJ.log("processing {0}".format(imageLabels['snapName']))
	images['nuclear'] = ImagePlus('nuclear', imp.getImageStack().getProcessor(parameters["nuclear"])).duplicate()
	images['bait'] = ImagePlus('bait', imp.getImageStack().getProcessor(parameters["bait"])).duplicate()
	"""use the DAPI channel to call nuclei"""
	DAPIoverlay = findnucleus(images['nuclear'], parameters)
	if not DAPIoverlay:
		IJ.log("no nuclei called in {0}".format(imageLabels['snapName']))
		return None
	totalnuclei = Overlay.size(DAPIoverlay)
	IJ.log("{0} nuclei found in {1}".format(totalnuclei, imageLabels['snapName']))
	baitoverlay = findarray(images, DAPIoverlay, totalnuclei, parameters)
	if not baitoverlay:
		IJ.log("no arrays coincident with called nuclei in {0}".format(imageLabels['snapName']))
		return None
	totalarray = Overlay.size(baitoverlay)
	IJ.log("{0} array(s) found in {1}".format(totalarray, imageLabels['snapName']))
	DAPIoverlay, totalnuclei = nucFilter(DAPIoverlay, baitoverlay, totalarray, totalnuclei)
	finalOverlay = nucArraypairer(DAPIoverlay,baitoverlay, totalarray, totalnuclei)
	if Overlay.size(finalOverlay) == 0:
		IJ.log("no single coincident arrays and nuclei identified in snap {0}".format(imageLabels['snapName']))
		return None
	table = measureImage(imp, finalOverlay, parameters, imageLabels)
	if len(table) > 0:
		"""save rois to output directory so can check success of array/nucleus caller and see which specific arrays & nuclei were identified"""
		roiSaver(finalOverlay, outputDir, imageLabels['snapName'])
		return table
	else:
		IJ.log("table not generated for {0} even though transfected cells overlapped specific nuclei.".format(imageLabels['snapName']))
		return None

def measureImage(imp, overlay, parameters, imageLabels):
	overlaySize = Overlay.size(overlay)
	imp.setOverlay(overlay)
	table = []
	channels = {"prey":parameters["prey"], "bait": parameters["bait"]}
	for label, channel in channels.items():
		imp.setSlice(channel)
		temp = []
		for i in range(overlaySize):
			roi = overlay.get(i)
			roi.setImage(imp)
			roi.setPosition(imp)
			Rname = Roi.getName(roi)
			roiStat = roi.getStatistics()
			temp.append(imageLabels['imagefile'])
			temp.append(imageLabels['date'])
			temp.append(imageLabels['snapName'])
			temp.append(label)
			temp.append(Rname)
			temp.append(roiStat.area)
			temp.append(roiStat.mean)
			temp.append(roiStat.median)
		table.append(temp)
	return table

def nucArraypairer(DAPIoverlay, baitoverlay, totalarray, totalnuclei):
	finalOverlay = Overlay()
	j = list(range(totalarray))
	for i in range(totalnuclei):
		tempNo = str(i + 1)
		roi = DAPIoverlay.get(i)
		nucPoints = roi.getContainedPoints()
		breakVal = 0
		for k, element in enumerate(j):
			roi2 = baitoverlay.get(element)
			arrayPoints = roi2.getContainedPoints()
			overlapTest = bool(set(arrayPoints) & set(nucPoints))
			if overlapTest == False:
				continue
			nucleoplasm = ShapeRoi(roi).not(ShapeRoi(roi2))
			Roi.setName(roi, "nucleus_" + tempNo)
			Roi.setName(roi2, "array_" + tempNo)
			Roi.setName(nucleoplasm, "nucleoplasm_" + tempNo)
			finalOverlay.add(roi)
			finalOverlay.add(roi2)
			finalOverlay.add(nucleoplasm)
			j.remove(element)
			break
	return finalOverlay

def nucFilter(DAPIoverlay, baitoverlay, totalarray, totalnuclei):
	"""remove all nuclei that do not contain precisely one array from DAPIoverlay"""
	completed = []
	j = list(range(totalarray))
	i = 0
	while i < totalnuclei:
		roi = DAPIoverlay.get(i)
		nucPoints = roi.getContainedPoints()
		overlapCounter = 0
		for k, element in enumerate(j):
			roi2 = baitoverlay.get(element)
			arrayPoints = roi2.getContainedPoints()
			overlapTest = bool(set(arrayPoints) & set(nucPoints))
			if overlapTest == True:
				j.remove(element)
				overlapCounter += 1
		if overlapCounter != 1:
			DAPIoverlay.remove(roi)
			totalnuclei -= 1
		else:
			i += 1
	return DAPIoverlay, totalnuclei

def processDirectory(inputDir):
	sep = System.getProperty("file.separator")
	outputDir = inputDir[0].getParent() + sep + "output" + sep
	cziFinder = re.compile(".*czi$")
	pathList = [image.getAbsolutePath() for directory in inputDir for image in directory.listFiles() if cziFinder.match(image.getAbsolutePath())]
	if pathList == []:
		IJ.log("No .czi files found in input directories. Exiting")
		IJ.error("No .czi files found in input directories. Exiting")
		return
	pathLen = len(pathList)
	CZIinfo = getCZIinfo(pathList[0])
	print CZIinfo
	if CZIinfo['SizeC'] < 2:
		IJ.log("F2H processing requires a minimum of 2 channels. Exiting")
		IJ.error("F2H processing requires a minimum of 2 channels. Exiting")
		return
	else:
		parameters = frameMaker().dialogBuilder(CZIinfo)
		if parameters == {}:
			return
		for key, value in parameters.items():
			print "{0}: {1}".format(key, value)
	if File(outputDir).exists() == False:
		File(outputDir).mkdir()
		IJ.log("Output directory created at {0}".format(outputDir))
	else:
		IJ.log("Output directory exists at {0}".format(outputDir))
	outputArray = []
	for i, image in enumerate(pathList):
		IJ.log("processing {0}. Image {1} out of {2}".format(image, i + 1, pathLen))
		outputArray.append(imageProcessor(image, parameters, outputDir))
		IJ.log("finished processing {0}. Image {1} out of {2}".format(image, i + 1, pathLen))
	IJ.log("image processing finished")
	log = IJ.getLog()
	resultsSaver(log, outputDir, "F2H_log", ".txt")
	table = resultsTablemaker(outputArray)
	if ResultsTable.size(table) == 0:
		IJ.log("No transfected cells found. Bye.")
		return
	resultsSaver(table, outputDir, "F2H_results", ".csv")
	WindowManager.closeAllWindows()

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
	elif extension == ".zip":
		rm = RoiManager.getInstance()
		rm.runCommand("save selected", filepath)
	IJ.log("Output saved to {0}".format(filepath))

def resultsTablemaker(outputArray):
	"""fill results table"""
	IJ.log("filling results table")
	table = ResultsTable()
	colnames = collections.deque(["Path", "Date", "Name", "Channel", "ROI", "Area", "Mean", "Median"])
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

def roiSaver(overlay, output, name):
	"""move overlay to roimanager and save"""
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.reset()
	for roi in overlay:
		rm.addRoi(roi)
	rm.deselect()
	resultsSaver("roi", output, name + "_rois", ".zip")
	#rm.runCommand("save selected", outputDir + name + "_rois.zip")

if __name__ in ["__builtin__", "__main__"]:
	processDirectory(inputDir)
