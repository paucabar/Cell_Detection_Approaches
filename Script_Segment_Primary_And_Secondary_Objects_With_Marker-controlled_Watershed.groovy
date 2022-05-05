#@ ImagePlus imp
#@ UpdateService updateService
#@ CommandService command
#@ ConvertService convertService

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
def impCyt = dup.run(imp, 1, 1, 1, 1, 1, 1);
def impNuc = dup.run(imp, 2, 2, 1, 1, 1, 1);

// run StarDist
def impNucLabels = runStarDist(impNuc, 0.5, 0.25)
setDisplayMinAndMax(impNucLabels)
impNucLabels.show()

// analyze regions
def ar = new AnalyzeRegions()
ar.setup("Area", impNucLabels)
def table = ar.process(impNucLabels)
table.show("Results")
int count = table.size()
println "$count labels"

// get values from the result table and store the index of labels that fail to meet one or more criteria
def labelDiscard = []
for (i in 0..count-1) {
	int area = table.getValue("Area", i)
	if (area < 105 || area > 2500) {
		labelDiscard.add(i+1)
	}
}
def labelDiscardFloat = labelDiscard as int[]
println "Discard labels $labelDiscardFloat"

// replace selected labels by 0
int replaceBy = 0
def ipImpNucLabels = impNucLabels.getProcessor()
rlv = new ReplaceLabelValues()
rlv.process(ipImpNucLabels, labelDiscardFloat, replaceBy)
//def impFilteredLabels = convertService.convert(ipImpNucLabels, ImagePlus.class)
def impFilteredLabels = new ImagePlus("Nuc Filtered Labels", ipImpNucLabels)
impFilteredLabels.show()

// marker-controlled watershed
impCytLabels = runMarkerControlledWatershed(impCyt, impNucLabels)
setDisplayMinAndMax(impCytLabels)
impCytLabels.show()

// replace -1.0 by 0.0
def replaceByFloat = replaceBy as float
def arrayMinusOne = [-1.0] as float[]
def ipImpCytLabels = impCytLabels.getProcessor()
rlv.process(ipImpCytLabels, arrayMinusOne, replaceByFloat)
def result = new ImagePlus("Result", ipImpCytLabels)
result.show()

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
	ic.convertToGray8()
	def gb = new GaussianBlur()
	def ipInputGaussian = inputGausian.getProcessor()

    gb.blurGaussian(ipInputGaussian, 1.5)
    ipInputGaussian.setAutoThreshold("Triangle dark")
	def ipBinaryMask = ipInputGaussian.createMask()
	mcwt = new MarkerControlledWatershedTransform2D (ipInput, ipLabels, ipBinaryMask, 8, 5.0)
	// I guess this connectivity setup  [int 8] establishes
	// diagonal connectivity (8-connected vs 4-connected)...
	
	labelsCell = mcwt.applyWithPriorityQueue() // label -1 for not detected cytoplasm
	def impCytLabels = new ImagePlus("Cyt Labels", labelsCell)
	return impCytLabels
}