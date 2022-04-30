#@ ImagePlus imp

impNuc = "nuc"
impCyt = "cyt"
impCytMask = "cytMask"
impNucLabels = "nucLabels"
impCytLabels = "cytLabels"

run("Duplicate...", "title="+impCyt+" duplicate channels=1");
selectImage(imp);
run("Duplicate...", "title="+impNuc+" duplicate channels=2");

selectWindow(impNuc);
run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+impNuc+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.5', 'nmsThresh':'0.25', 'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
selectWindow("Label Image");
rename(impNucLabels);
filterLabels(impNucLabels, 105, 2500);

selectWindow(impCyt);
run("Duplicate...", "title="+impCytMask);
run("Gaussian Blur...", "sigma=1.50");
run("8-bit");
setAutoThreshold("Triangle dark");
run("Convert to Mask");
run("Marker-controlled Watershed", "input="+impCyt
	+" marker="+impNucLabels+" mask="+impCytMask
	+" compactness=5 calculate use");
rename(impCytLabels);

selectWindow(impCytLabels);
run("Label image to ROIs", "rm=[RoiManager[size=102, visible=true]]");
roiManager("Show None");
selectImage(imp);
roiManager("Show All without labels");

close(impNuc);
close(impCyt);
close(impCytMask);
selectWindow("Log");
run("Close");
run("Tile");


function filterLabels(imageTitle, minArea, maxArea) {
	selectWindow(imageTitle);
	run("Analyze Regions", "area");
	IJ.renameResults(impNucLabels+"-Morphometry", "Results");
	setOption("ExpandableArrays", true);
	labelDiscard=newArray();
	nLabels=nResults;
	discardCount=0;
	for (i=0; i<nLabels; i++) {
		area=getResult("Area", i);
		if (area < minArea || area > maxArea) {
			labelDiscard[discardCount]=i+1;
			discardCount++;
		}
	}
	if (discardCount != 0) {
		discardString="";
		for (i=0; i<labelDiscard.length; i++) {
			n=d2s(labelDiscard[i], 0);
			discardString+=n;
			if(i<labelDiscard.length-1) {
				discardString+=",";
			}
		}
		print(discardString);
		run("Replace/Remove Label(s)", "label(s)="+discardString+" final=0");
	}
	selectWindow("Results");
	run("Close");
}
