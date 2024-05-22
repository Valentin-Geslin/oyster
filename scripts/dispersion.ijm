path = File.getDefaultDir;
root_path = path+"images_histoqc_masked"+File.separator;
dispersion_path = path+"results"+File.separator+"dispersion"+File.separator;

function dispesion(imagename, distance, known) { 
	close("*");
	file = root_path+imagename+".png";
	
	if (File.exists(file)) {
		open(file);
		run("Set Scale...", "distance="+distance+" known="+known+" unit=µm");
		run("16-bit");
		setAutoThreshold("Default dark");
		run("Threshold...");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Particle Distribution (2D)", "min_size=0.1 max_size=Infinity min_circularity=0 max_circularity=1 neighbor=[centroid NND] statistical=mean conficence=95%");
		if (!File.isDirectory(dispersion_path)) {
			File.makeDirectory(dispersion_path);
		}
		saveAs("Text", dispersion_path+imagename+".txt");
		run("Close");
		close();
	}
}

dispesion("64851", "3720", "27017.77");
dispesion("64852", "4320", "27071.77");
dispesion("64853", "3840", "27939.97");
dispesion("64854", "3840", "27939.97");
dispesion("64856", "3720", "27066.85");
dispesion("64857", "3960", "28813.10");
dispesion("16049-2-3101002", "3600", "26193.72");
dispesion("16049-3-3201002", "3720", "27066.85");
dispesion("16049-5-16110141B00201002", "3480", "25325.21");
dispesion("16049-6-3301002", "3720", "27066.85");
dispesion("16049-7-16110141B00301002", "3720", "27071.77");
dispesion("16049-8-3401002", "3480", "25320.60");
dispesion("16049-10-16110141B00401002", "3720", "27071.77");
dispesion("16045-1-101002", "4080", "29686.22");
dispesion("16045-2-16110141B00901002", "3840", "27945.06");
dispesion("16045-3-16110141B01001002", "3840", "27945.06");
dispesion("16045-4-201002", "3840", "27939.97");
dispesion("16045-5-16110141B01101002", "4440", "32311.47");
dispesion("16045-6-301002", "3600", "26193.72");
dispesion("16045-7-401003", "4320", "31432.47");
dispesion("16045-8-16110141B01201002", "3960", "28818.34");
dispesion("16045-9-16110141B01301002", "3312", "24091.65");
dispesion("16045-10-16110141B01401002", "3480", "25325.21");
dispesion("16045-11-16110141B01501002", "3480", "25325.21");
dispesion("16045-13-601002", "2720", "19780.02");
dispesion("16045-14-16110141B01601002", "3720", "27071.77");
dispesion("16045-7-701002", "4080", "29686.22");
dispesion("16045-17-801002", "3960", "28813.10");
dispesion("16045-18-901002", "4200", "30559.35");
dispesion("16045-19-1001002", "3960", "28813.10");
dispesion("16045-20-1101002", "3360", "24447.48");
dispesion("16045-21-1201002", "3240", "23574.35");
dispesion("16045-22-16110141B01701002", "2520", "18338.94");
dispesion("16045-23-16110141B01801002", "3240", "23578.64");
dispesion("16045-24-16110141B01901002", "3480", "25325.21");
dispesion("16045-25-1301002", "3480", "25320.60");
dispesion("16045-28-16110141B02001002", "3008", "21880.34");
dispesion("16045-29-16110141B02101002", "2880", "20958.79");
dispesion("16045-31-16110141B02201002", "3960", "28818.34");
dispesion("16045-32-1701002", "4080", "29686.22");
dispesion("16045-50-16110141B02401002", "3360", "24451.92");