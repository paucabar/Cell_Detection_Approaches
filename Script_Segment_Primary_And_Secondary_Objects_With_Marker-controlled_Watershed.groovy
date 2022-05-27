#@ ImagePlus imp
#@ UpdateService updateService
#@ UIService ui
#@ CommandService command
#@ ConvertService convertService
#@ String (label=" ", value="Channels", visibility=MESSAGE, persist=false) message1
#@ Integer (label="Nuclei Marker", value=2, max=4, min=1, style="listBox") nucleiChannel
#@ Integer (label="Cell Marker", value=1, max=4, min=1, style="listBox") cellChannel
#@ String (label=" ", value="StarDist", visibility=MESSAGE, persist=false) message2
#@ Double (label="StarDist Score Threshold", value=0.5, max=1.0, min=0.0, stepSize=0.05, style="slider") scoreSD
#@ Double (label="StarDist Overlap Threshold", value=0.25, max=1.00, min=0.00, stepSize=0.05, style="slider") overlapSD
#@ String (label=" ", value="Nuclei Size", visibility=MESSAGE, persist=false) message3
#@ Integer (label="Min Area", value=105, style="listBox") minArea
#@ Integer (label="Max Area", value=2500, style="listBox") maxArea
#@ String (label=" ", value="Cell Marker Mask", visibility=MESSAGE, persist=false) message4
#@ Double (label="Enhance Contrast [% Sat Pixels]", value=3.0, stepSize=0.5, style="listBox") saturatedPixels
#@ Double (label="Gaussian Blur [radius]", value=1.5, stepSize=0.5, style="listBox") gaussianRadius
#@ String (label="Thresholding Method", choices={"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, value="Triangle", style="listBox") thresholdingMethod
#@ String (label=" ", value="Marker-controlled Watershed", visibility=MESSAGE, persist=false) message5
#@ String (label="Connectivity", choices={"4", "8"}, value=8, style="radioButtonHorizontal") connectivityString
#@ Double (label="Compactness", value=5.0, stepSize=0.5, style="listBox") compactness

import ij.IJ
import ij.plugin.Duplicator
import de.csbdresden.stardist.StarDist2D
import ij.ImagePlus
import inra.ijpb.plugins.AnalyzeRegions
import inra.ijpb.label.edit.ReplaceLabelValues
import ij.plugin.filter.GaussianBlur
import inra.ijpb.watershed.MarkerControlledWatershedTransform2D
import ij.process.ImageConverter

// check update sites
boolean checkStarDist = isUpdateSiteActive("StarDist");
boolean checkCSBDeep = isUpdateSiteActive("CSBDeep");
boolean checkMorphoLibJ = isUpdateSiteActive("IJPB-plugins");

// exit if any update site is missing
if (!checkStarDist || !checkCSBDeep || !checkMorphoLibJ) {
	return
}

// duplicate channels
def dup = new Duplicator()
def impCyt = dup.run(imp, cellChannel, cellChannel, 1, 1, 1, 1);
def impCyt2 = dup.run(imp, cellChannel, cellChannel, 1, 1, 1, 1);
def impNuc = dup.run(imp, nucleiChannel, nucleiChannel, 1, 1, 1, 1);

// run StarDist
def impNucLabels = runStarDist(impNuc, scoreSD, overlapSD)
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
def impFilteredLabels = new ImagePlus("Nuclei Labels", ipImpNucLabels)
impFilteredLabels.show()

// cytoplasm semantic segmentation
impCytoplasmMask = semanticSegmentation(impCyt2)
//impCytoplasmMask.show()

// marker-controlled watershed
def connectivity = connectivityString as int
println "Connectivity: $connectivity-connected"
impCytLabels = runMarkerControlledWatershed(impCyt, impFilteredLabels, impCytoplasmMask, connectivity)
setDisplayMinAndMax(impCytLabels)
//impCytLabels.show()

// in cell labels: replace -1.0 by 0.0
def replaceByFloat = replaceBy as float
def arrayMinusOne = [-1.0] as float[]
def ipImpCytLabels = impCytLabels.getProcessor()
rlv.process(ipImpCytLabels, arrayMinusOne, replaceByFloat)
def impCellLabels = new ImagePlus("Cytoplasm Labels", ipImpCytLabels)
def ic = new ImageConverter(impCellLabels)
ic.setDoScaling(false)
ic.convertToGray16()
ic.setDoScaling(true)
setDisplayMinAndMax(impCellLabels)
impCellLabels.show()

def isUpdateSiteActive (updateSite) {
	checkUpdate = true
	if (! updateService.getUpdateSite(updateSite).isActive()) {
    	ui.showDialog "Please activate the $updateSite update site"
    	checkUpdate = false
	}
	return checkUpdate
}

def runStarDist(input, scoreTh, overlapTh) {
	res = command.run(StarDist2D, false,
			"input", input,
			"modelChoice", "Versatile (fluorescent nuclei)",
			"probThresh", scoreTh,
			"nmsThresh", overlapTh,
			"outputType", "Label Image",
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

def semanticSegmentation(input) {
	IJ.run(input, "Enhance Contrast...", "saturated=$saturatedPixels update")
	IJ.run(input, "Apply LUT", "");
	def gb = new GaussianBlur()
	def ipInput = input.getProcessor()
	gb.blurGaussian(ipInput, gaussianRadius)
	ipInput.setAutoThreshold("$thresholdingMethod dark")
	def ipBinaryMask = ipInput.createMask()
	def impBinaryMask = new ImagePlus("Cyt Mask", ipBinaryMask)
	return impBinaryMask
}

def runMarkerControlledWatershed(input, labels, mask, con_4or8) {
	def ipInput = input.getProcessor()
	def ipLabels = labels.getProcessor()
	def ipMask = mask.getProcessor()
	mcwt = new MarkerControlledWatershedTransform2D (ipInput, ipLabels, ipMask, con_4or8, compactness)
	labelsCell = mcwt.applyWithPriorityQueue() // label -1 for not detected cytoplasm
	def impCytLabels = new ImagePlus("Cyt Labels", labelsCell)
	return impCytLabels
}