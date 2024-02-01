#@ UpdateService updateService
#@ UIService ui
#@ CommandService command
#@ ConvertService convertService
#@ File(label="Image File", style="open") file
#@ File(label="StarDist Model", style="open") modelFile
#@ String (label=" ", value="Channels", visibility=MESSAGE, persist=false) message1
#@ Integer (label="Nuclear Marker", value=1, max=4, min=1, style="listBox") nucleiChannel
#@ Integer (label="Membrane Marker", value=2, max=4, min=1, style="listBox") membraneChannel
#@ Integer (label="Slice", value=1, max=10, min=1, style="listBox") slice
#@ String (label=" ", value="StarDist", visibility=MESSAGE, persist=false) message2
#@ Double (label="StarDist Score Threshold", value=0.5, max=1.0, min=0.0, stepSize=0.05, style="slider") scoreSD
#@ Double (label="StarDist Overlap Threshold", value=0.25, max=1.00, min=0.00, stepSize=0.05, style="slider") overlapSD
#@ String (label=" ", value="Nuclei Size", visibility=MESSAGE, persist=false) message3
#@ Integer (label="Min Area", value=105, style="listBox") minArea
#@ Integer (label="Max Area", value=2500, style="listBox") maxArea
#@ Integer (label="Erode Iterations", value=1, max=99, min=1, style="listBox") iterations
#@ String (label="LUT", choices={"mpl-viridis", "mpl-plasma", "mpl-inferno", "mpl-magma"}, value="mpl-viridis", style="listBox") lutName


import ij.io.Opener
import ij.IJ
import ij.process.LUT
import java.awt.Color
import ij.plugin.Duplicator
import de.csbdresden.stardist.StarDist2D
import ij.ImagePlus
import inra.ijpb.plugins.AnalyzeRegions
import inra.ijpb.label.edit.ReplaceLabelValues
import inra.ijpb.binary.distmap.ChamferMask2D
import inra.ijpb.label.filter.ChamferLabelErosion2DShort
import inra.ijpb.watershed.MarkerControlledWatershedTransform2D
import ij.plugin.filter.GaussianBlur
import ij.plugin.ImageCalculator
import ij.plugin.ContrastEnhancer
import ij.process.ImageConverter
import ij.process.ImageProcessor
import ij.process.ByteProcessor
import ij.Prefs
import ij.plugin.frame.RoiManager
import ij.gui.Roi
import ij.plugin.filter.ThresholdToSelection
import ij.measure.ResultsTable
import ij.IJ
import ij.plugin.LutLoader 
import ij.process.LUT


def isUpdateSiteActive (updateSite) {
	checkUpdate = true
	if (! updateService.getUpdateSite(updateSite).isActive()) {
    	ui.showDialog "Please activate the $updateSite update site"
    	checkUpdate = false
	}
	return checkUpdate
}

ImagePlus openImage(File file) {
	String path = file.getAbsolutePath()
	def opener = new Opener()
	String extension = path[path.lastIndexOf('.')+1..-1]
	println "Importing $extension file"
	ImagePlus imp = opener.openUsingBioFormats(path)


	return imp
}

// set composite LUTs
void setLUTs(ImagePlus imp) {
	colorList = [Color.RED, Color.GREEN, Color.BLUE, Color.MAGENTA, Color.CYAN, Color.YELLOW]
	luts = imp.getLuts()
	int nChannels = imp.nChannels
	for (i in 0..nChannels-1) {
		luts[i] = LUT.createLutFromColor(colorList[i])
	}
	
	imp.setLuts(luts)
	int displayMode = imp.getDisplayMode()
	if (displayMode != IJ.COMPOSITE) {
		imp.setDisplayMode(IJ.COMPOSITE)
	}
}

String getModelPath(File modelFile) {
	// Get the file path as a string
	String filePath = modelFile.getAbsolutePath()
	println("Original File Path: $filePath")
	
	// Replace file separators with four file separators
	String modifiedPath = filePath.replace(File.separator, File.separator * 4)
	println("Modified File Path: $modifiedPath")
	return modifiedPath
}

ImagePlus runStarDist(input, scoreTh, overlapTh, modelFile) {
	res = command.run(StarDist2D, false,
			"input", input,
			"modelChoice", "Model (.zip) from File",
			"probThresh", scoreTh,
			"nmsThresh", overlapTh,
			"outputType", "Label Image",
			"modelFile", getModelPath(modelFile),
			).get()
	def inputLabel = res.getOutput("label")
	def impInputLabel = convertService.convert(inputLabel, ImagePlus.class)
	return impInputLabel
}

void setDisplayMinAndMax(image) {
	def ip = image.getProcessor()
	maxImage = ip.getStats().max as int
	println "Set display 0 - $maxImage"	
	image.setDisplayRange(0, maxImage)
	IJ.run(image, "glasbey inverted", "")
}

//ImagePlus labelToBinary(ImagePlus imp) {
//	ImageProcessor ip = imp.getProcessor()
//	ip.setThreshold (1, 255, ImageProcessor.NO_LUT_UPDATE)
//	ImageProcessor ipBinary = ip.createMask() // image processor
//	ImagePlus impBinary = new ImagePlus("Binary Mask", ipBinary)
//	return impBinary
//}

ImagePlus erodeLabels(ImagePlus imp, double radius) {
	def labelErosion = new ChamferLabelErosion2DShort(ChamferMask2D.CHESSBOARD, radius)
	ImageProcessor ipEroded = labelErosion.process(imp.getProcessor())
	ImagePlus impEroded = new ImagePlus("Label Eroded", ipEroded)
	return impEroded
}

ImagePlus subtractLabels(ImagePlus imp1, ImagePlus imp2) {
	ImageCalculator imgCalculator = new ImageCalculator()
	ImagePlus result = imgCalculator.run(imp1, imp2, "subtract")
	return result
}

// get rois from labels
RoiManager labelsToRois(ImagePlus imp, RoiManager rm) {
	ImageProcessor ip = imp.getProcessor()
	ThresholdToSelection getSelection = new ThresholdToSelection()
	def max = imp.getStatistics().max
	for (i in 1..max) {
		ip.setThreshold(i, i, ImageProcessor.NO_LUT_UPDATE)
	    Roi roiTemp = getSelection.convert(ip)
	    if (roiTemp != null) {
			//roiTemp.setColor(Color.RED)
			roiTemp.setName(String.format("%03d", i))
			roiTemp.setStrokeWidth(1)
			rm.addRoi(roiTemp)
	    }
	}
	return rm
}

// This function calculates the specified number of bins (each represented as a tuple with minimum and maximum values) based on a list of numeric values
Collection getBins(List<Number> values, int numBins) {
	// Calculate bin width
	def binWidth = (values.max() - values.min()) / numBins
	
	// Create bins as a list of tuples
	def bins = (0..<numBins).collect { index ->
	    def min = values.min() + index * binWidth
	    def max = index == numBins - 1 ? values.max() + 1e-6 : min + binWidth
	    [min, max]
	}
	return bins
}

// Function to find the bin for a new value and return a tuple representing the bin's min and max values
// Returns: Tuple
def findBinForValue(value, bins) {
    return bins.find { bin ->
        value >= bin[0] && value < bin[1]
    }
}

// Function to label values with bins
List<Integer> labeledValues(List<Number> values, int numBins) {
    def bins = getBins(values, numBins)
    
    // Create a list to store the labels
    def labels = []

    // Iterate through the values and label them with the corresponding bin
    values.each { value ->
        def bin = findBinForValue(value, bins)
        if (bin) {
            // Add 1 to the bin index to start from 1
            labels.add(bins.indexOf(bin) + 1)
        } else {
            // Handle values that don't belong to any bin
            labels.add(0) // You can use 0 or another value to represent this case
        }
    }
    return labels
}

def getLutFromName(String name) {
	// Get the ImageJ directory
	def imagejDirectory = IJ.getDirectory("imagej")
	
	// Specify the LUT file path relative to the ImageJ directory
	String lutDirectory = "luts"
	String lutFileName = name + ".lut"
	String pathLut = imagejDirectory + lutDirectory + File.separator + lutFileName
	
	// Load LUT
	LUT originalLut = LutLoader.openLut(pathLut)
	
	// Split the list into 3 lists with 256 values each
	bytesList = originalLut.getBytes()
	def splitLists = []
	def batchSize = 256
	for (int i = 0; i < bytesList.size(); i += batchSize) {
	    def endIndex = Math.min(i + batchSize, bytesList.size())
	    def sublist = bytesList[i..endIndex - 1] as byte[]
	    splitLists << sublist
	}
	
	// Replace the first value of each splitList with 0 as a byte
	for (def i = 0; i < splitLists.size(); i++) {
	    splitLists[i][0] = (byte) 0
	}
	
	// Use modified bytes to create new LUT
	LUT newLUT = new LUT(splitLists[0], splitLists[1], splitLists[2])
	return newLUT
}

/**
 * Creates a label image from the RoiManager
 * Uses Roi indexes as labels
 */
ImagePlus resultColormap(ImagePlus imp, RoiManager rm, List<Integer> binLabelsList, String lutName) {
    impColormap = IJ.createImage("Labeling", "32-bit black", imp.getWidth(), imp.getHeight(), 1)
    ip = impColormap.getProcessor()
    rm.getRoisAsArray().eachWithIndex { roi, index ->
        ip.setColor(binLabelsList[index]) // set membrane intensity as color
        ip.fill(roi)
    }
    ip.resetMinAndMax()
    LUT lut = getLutFromName(lutName)
    impColormap.setLut(lut)
    return impColormap
}

// check update sites
boolean checkStarDist = isUpdateSiteActive("StarDist");
boolean checkCSBDeep = isUpdateSiteActive("CSBDeep");
boolean checkMorphoLibJ = isUpdateSiteActive("IJPB-plugins");

// exit if any update site is missing
if (!checkStarDist || !checkCSBDeep || !checkMorphoLibJ) {
	return
}

//import image
imp = openImage(file)
imp.show()
setLUTs(imp)

// some checks on imp
dims = imp.getDimensions() // default order: XYCZT
if (dims[2] < nucleiChannel) {
	ui.showDialog "The nuclear channel has been set to $nucleiChannel, while the image contains only ${dims[2]} channels."
	return
} else if (dims[2] < membraneChannel) {
	ui.showDialog "The membrane channel has been set to $membraneChannel, while the image contains only ${dims[2]} channels."
	return
} else if (nucleiChannel == membraneChannel) {
	ui.showDialog "WARNING: Both the membrane channel and the nuclear channel have been set to the same value: [$membraneChannel]."
} else if (dims[3] < slice) {
	ui.showDialog "The slice has been set to $slice, while the image contains only ${dims[3]} slices."
	return
}

// duplicate channels
Duplicator duplicator = new Duplicator()
ImagePlus impMembrane = duplicator.run(imp, membraneChannel, membraneChannel, slice, slice, 1, 1)
ImagePlus impNuc = duplicator.run(imp, nucleiChannel, nucleiChannel, slice, slice, 1, 1)
//impMembrane.show()
//impNuc.show()

// run StarDist
def impNucLabels = runStarDist(impNuc, scoreSD, overlapSD, modelFile)
setDisplayMinAndMax(impNucLabels)
//impNucLabels.show()

// analyze regions
def ar = new AnalyzeRegions()
//ar.setup("Area", impNucLabels)
def table = ar.process(impNucLabels)
//table.show("Results")
int count = table.size()
println "$count labels"

// get values from the result table and store the index of labels
// that fail to meet one or more criteria
def labelDiscard = []
for (i in 0..count-1) {
	int area = table.getValue("Area", i)
	if (area < minArea || area > maxArea) {
		labelDiscard.add(i+1)
	}
}
def labelDiscardInt = labelDiscard as int[]
println "Discard labels $labelDiscardInt"

// in nuclei labels: replace selected labels by 0
int replaceBy = 0
def ipImpNucLabels = impNucLabels.getProcessor()
rlv = new ReplaceLabelValues()
rlv.process(ipImpNucLabels, labelDiscardInt, replaceBy)
ImagePlus impFilteredLabels = new ImagePlus("Nuclei Labels", ipImpNucLabels)
//impFilteredLabels.show()

// get the nuclei membrane and inner labels
ImagePlus impEroded = erodeLabels(impFilteredLabels, 3.0)
ImagePlus impBorder = duplicator.run(impFilteredLabels, 1, 1, 1, 1, 1, 1)
subtractLabels(impBorder, impEroded)
//impEroded.show()
//impBorder.show()

// create RoiManagers
RoiManager nucleusRoiManager = new RoiManager(false)
RoiManager membraneRoiManager = new RoiManager(false)
RoiManager innerRoiManager = new RoiManager(false)

// fill RoiManagers
nucleusRoiManager = labelsToRois(impFilteredLabels, nucleusRoiManager)
membraneRoiManager = labelsToRois(impBorder, membraneRoiManager)
innerRoiManager = labelsToRois(impEroded, innerRoiManager)
assert nucleusRoiManager.getCount() == membraneRoiManager.getCount() && membraneRoiManager.getCount() == innerRoiManager.getCount(), "The count of the nuclei, nuclear membrane, and nuclear intern ROIs doesn't match. Perhaps you eroded too much??"

// measure
resultsMap = [ : ]
for (i in 0..nucleusRoiManager.getCount()-1) {
	// get the 3 Rois and check they correspond to the same label
	nucRoi = nucleusRoiManager.getRoi(i)
	memRoi = membraneRoiManager.getRoi(i)
	inRoi = innerRoiManager.getRoi(i)
	assert nucRoi.getName() == memRoi.getName() && memRoi.getName() == inRoi.getName(), "Roi names don't match"
	// measure mean intensity on membrane staining channel 
	impMembrane.setRoi(memRoi)
	memMean = memRoi.getStatistics().mean
	impMembrane.setRoi(inRoi)
	inMean = inRoi.getStatistics().mean
	resultsMap[nucRoi.getName()] = [memMean, inMean, memMean/inMean]
}

// Get the keys and values from the map
def keysList = resultsMap.keySet()
def memMeanList = resultsMap.collect { key, value -> value[0] }
def inMeanList = resultsMap.collect { key, value -> value[1] }
def ratioList = resultsMap.collect { key, value -> value[2] }


// create and fill result table
ResultsTable rt = new ResultsTable(keysList.size())
rt.setPrecision(5)
rt.setValues("Membrane Mean", memMeanList as double[])
rt.setValues("Internal Mean", inMeanList as double[])
rt.setValues("Mem/Int Rate", ratioList as double[])
// set table labels
for (i in 0..keysList.size()-1) {
    rt.setLabel(keysList[i], i)
}

// show results table
rt.show("Results Table")

// create colormap
int nBins = 20
Collection bins = getBins(memMeanList, nBins)
List<Integer> binLabelsList = labeledValues(memMeanList, nBins)
colormap = resultColormap(impFilteredLabels, nucleusRoiManager, binLabelsList, lutName)
colormap.show()
impMembrane.show()
IJ.run(impMembrane, "Grays", "")
IJ.run("Add Image...", "image=" + colormap.getTitle() + " x=0 y=0 opacity=35 zero")
IJ.run("Select None")
return