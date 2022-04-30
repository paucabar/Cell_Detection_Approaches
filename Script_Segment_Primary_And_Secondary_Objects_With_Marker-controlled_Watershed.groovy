#@ ImagePlus imp
#@ UpdateService updateService
#@ CommandService command
#@ ConvertService convertService

import ij.IJ
import ij.plugin.Duplicator
import de.csbdresden.stardist.StarDist2D
import ij.ImagePlus
import inra.ijpb.plugins.AnalyzeRegions
import inra.ijpb.watershed.MarkerControlledWatershedTransform2D

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
def impLabels = runStarDist(impNuc, 0.5, 0.25)
setDisplayMinAndMax(impLabels)
impLabels.show()

// analyze regions
ar = new AnalyzeRegions()
//def setup = ar.setup("area", impLabels)
def table = ar.process(impLabels)
table.show("Results")
int count = table.size()
println "$count labels"

// get values from the result table and store the index of labels that fail to meet one or more criteria
def labelDiscard = []
for (i in 0..count-1) {
	int area = table.getValue("Area", i)
	if (area < 120 || area > 2500) {
		labelDiscard.add(i+1)
	}
}
println "Discard labels $labelDiscard"

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