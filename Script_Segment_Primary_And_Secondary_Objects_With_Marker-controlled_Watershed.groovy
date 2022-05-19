#@ ImagePlus imp
#@ UpdateService updateService
#@ CommandService command
#@ ConvertService convertService
#@ Integer (label="Nuclei Channel", value=2, max=4, min=1, style="listBox") nucleiChannel
#@ Integer (label="Cell Channel", value=1, max=4, min=1, style="listBox") cellChannel
#@ Double (label="StarDist Score Threshold", value=0.5, max=1.0, min=0.0, stepSize=0.05, style="slider") scoreSD
#@ Double (label="StarDist Overlap Threshold", value=0.25, max=1.00, min=0.00, stepSize=0.05, style="slider") overlapSD
#@ Integer (label="Min Area", value=105, style="listBox") minArea
#@ Integer (label="Max Area", value=2500, style="listBox") maxArea
#@ Double (label="Gaussian Blur [radius]", value=1.5, stepSize=0.5, style="listBox") gaussianRadius
#@ String (label="Thresholding Method", choices={"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, value="Triangle", style="listBox") thresholdingMethod

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
boolean checkBIOP = isUpdateSiteActive("PTBIOP");

// exit if any update site is missing
// TODO

// duplicate channels
def dup = new Duplicator()
def impCyt = dup.run(imp, cellChannel, cellChannel, 1, 1, 1, 1);
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

// marker-controlled watershed
impCytLabels = runMarkerControlledWatershed(impCyt, impFilteredLabels)
setDisplayMinAndMax(impCytLabels)
//impCytLabels.show()

// in cell labels: replace -1.0 by 0.0
def replaceByFloat = replaceBy as float
def arrayMinusOne = [-1.0] as float[]
def ipImpCytLabels = impCytLabels.getProcessor()
rlv.process(ipImpCytLabels, arrayMinusOne, replaceByFloat)
def impCellLabels = new ImagePlus("Cytoplasm Labels", ipImpCytLabels)
def ic3 = new ImageConverter(impCellLabels)
ic3.setDoScaling(false)
ic3.convertToGray16()
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

def runMarkerControlledWatershed(input, labels) {
	def ipInput = input.getProcessor()
	def ipLabels = labels.getProcessor()
	
	def inputGausian = input.duplicate()
	def ic = new ImageConverter(inputGausian)
	ic.setDoScaling(true)
	ic.convertToGray8()
	def gb = new GaussianBlur()
	def ipInputGaussian = inputGausian.getProcessor()

    gb.blurGaussian(ipInputGaussian, gaussianRadius)
    ipInputGaussian.setAutoThreshold("$thresholdingMethod dark")
	def ipBinaryMask = ipInputGaussian.createMask()
	mcwt = new MarkerControlledWatershedTransform2D (ipInput, ipLabels, ipBinaryMask, 8, 5.0)
	// I guess this connectivity setup  [int 8] establishes
	// diagonal connectivity (8-connected vs 4-connected)...
	
	labelsCell = mcwt.applyWithPriorityQueue() // label -1 for not detected cytoplasm
	def impCytLabels = new ImagePlus("Cyt Labels", labelsCell)
	return impCytLabels
}