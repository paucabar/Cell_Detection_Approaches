#@ ImagePlus imp
#@ UpdateService updateService
#@ CommandService command
#@ ConvertService convertService

import ij.IJ
import ij.plugin.Duplicator
import de.csbdresden.stardist.StarDist2D
import ij.ImagePlus

// check update sites
a = isUpdateSiteActive("StarDist");
b = isUpdateSiteActive("CSBDeep");
c = isUpdateSiteActive("IJPB-plugins");
d = isUpdateSiteActive("PTBIOP");

// duplicate channels
def dup = new Duplicator()
def impCyt = dup.run(imp, 1, 1, 1, 1, 1, 1);
def impNuc = dup.run(imp, 2, 2, 1, 1, 1, 1);

// run StarDist
def impLabels = runStarDist(impNuc, 0.5, 0.25)
setDisplayMinAndMax(impLabels)
impLabels.show()

// analyze regions
IJ.run(impLabels, "Analyze Regions", "area");


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