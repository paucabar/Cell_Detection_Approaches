#@ ImagePlus imp
#@ UpdateService updateService
#@ UIService ui
#@ CommandService command
#@ ConvertService convertService
#@ File(label="StarDist Model", style="open") modelFile
#@ String (label=" ", value="Channels", visibility=MESSAGE, persist=false) message1
#@ Integer (label="Nuclei Marker", value=1, max=4, min=1, style="listBox") nucleiChannel
#@ Integer (label="Membrane Marker", value=2, max=4, min=1, style="listBox") membraneChannel
#@ String (label=" ", value="StarDist", visibility=MESSAGE, persist=false) message2
#@ Double (label="StarDist Score Threshold", value=0.5, max=1.0, min=0.0, stepSize=0.05, style="slider") scoreSD
#@ Double (label="StarDist Overlap Threshold", value=0.25, max=1.00, min=0.00, stepSize=0.05, style="slider") overlapSD
#@ String (label=" ", value="Nuclei Size", visibility=MESSAGE, persist=false) message3
#@ Integer (label="Min Area", value=105, style="listBox") minArea
#@ Integer (label="Max Area", value=2500, style="listBox") maxArea
#@ Integer (label="Erode Iterations", value=1, max=99, min=1, style="listBox") iterations

import ij.IJ
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


def isUpdateSiteActive (updateSite) {
	checkUpdate = true
	if (! updateService.getUpdateSite(updateSite).isActive()) {
    	ui.showDialog "Please activate the $updateSite update site"
    	checkUpdate = false
	}
	return checkUpdate
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

def runStarDist(input, scoreTh, overlapTh, modelFile) {
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

def setDisplayMinAndMax(image) {
	def ip = image.getProcessor()
	maxImage = ip.getStats().max as int
	println "Set display 0 - $maxImage"	
	image.setDisplayRange(0, maxImage)
	IJ.run(image, "glasbey inverted", "")
}

// Implements the Erode, Dilate, Open and Close commands in the Process/Binary submenu. 
public void run (ImageProcessor ip, String arg) {
    int fg = Prefs.blackBackground ? 255 : 0;
    foreground = ip.isInvertedLut() ? 255-fg : fg;
    background = 255 - foreground;
    ip.setSnapshotCopyMode(true);

	if (arg.equals("erode") || arg.equals("dilate")) {
        doIterations((ByteProcessor)ip, arg);
	} else if (arg.equals("open")) {
        doIterations(ip, "erode");
        doIterations(ip, "dilate");
    } else if (arg.equals("close")) {
        doIterations(ip, "dilate");
        doIterations(ip, "erode");
    }
    ip.setSnapshotCopyMode(false);
    ip.setBinaryThreshold();
}

void doIterations (ImageProcessor ip, String mode) {
    for (int i=0; i<iterations; i++) {
        if (Thread.currentThread().isInterrupted()) return;
        if (IJ.escapePressed()) {
            escapePressed = true;
            ip.reset();
            return;
        }
        if (mode.equals("erode"))
            ((ByteProcessor)ip).erode(count, background);
        else
            ((ByteProcessor)ip).dilate(count, background);
    }
}

ImagePlus labelToBinary(ImagePlus imp) {
	ImageProcessor ip = imp.getProcessor()
	ip.setThreshold (1, 255, ImageProcessor.NO_LUT_UPDATE)
	ImageProcessor ipBinary = ip.createMask() // image processor
	ImagePlus impBinary = new ImagePlus("Binary Mask", ipBinary)
	return impBinary
}

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
		//roiTemp.setColor(Color.RED)
		roiTemp.setName(String.format("%03d", i))
		roiTemp.setStrokeWidth(1)
		rm.addRoi(roiTemp)
	}
	return rm
}


// check update sites
boolean checkStarDist = isUpdateSiteActive("StarDist");
boolean checkCSBDeep = isUpdateSiteActive("CSBDeep");
boolean checkMorphoLibJ = isUpdateSiteActive("IJPB-plugins");

// exit if any update site is missing
if (!checkStarDist || !checkCSBDeep || !checkMorphoLibJ) {
	return
}

// duplicate channels
Duplicator duplicator = new Duplicator()
ImagePlus impMembrane = duplicator.run(imp, membraneChannel, membraneChannel, 1, 1, 1, 1);
ImagePlus impNuc = duplicator.run(imp, nucleiChannel, nucleiChannel, 1, 1, 1, 1);

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

/**
 * Creates a label image from the RoiManager
 * Uses Roi indexes as labels
 */
ImagePlus resultColormap(ImagePlus imp, RoiManager rm, map) {
    impColormap = IJ.createImage("Labeling", "32-bit black", imp.getWidth(), imp.getHeight(), 1)
    ip = impColormap.getProcessor()
    rm.getRoisAsArray().eachWithIndex { roi, index ->
        roiCode = roi.getName()
        metrics = map[roiCode]
        ip.setColor(metrics[0]) // set membrane intensity as color
        ip.fill(roi)
    }
    ip.resetMinAndMax()
    IJ.run(impColormap, "mpl-plasma", "")
    return impColormap
}

colormap = resultColormap(impFilteredLabels, nucleusRoiManager, resultsMap)
colormap.show()
impMembrane.show()
IJ.run(impMembrane, "Grays", "")
return