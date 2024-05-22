path = File.getDefaultDir;
root_path = path+"images_histoqc_masked"+File.separator;
mask_path = path+"mask_histoqc"+File.separator+"mask_use"+File.separator;

function colocalization(imagename, distance, known) {
	close("*");
	file = root_path+imagename+".png";
	mask = mask_path+imagename+".png";
	
	if (File.exists(file) && File.exists(mask)){
		open(file);
		rename(imagename+"_IHC_mask.png");
		run("Set Scale...", "distance="+distance+" known="+known+" unit=µm");
		run("16-bit");
		run("Threshold...");
		setOption("BlackBackground", false);
		setThreshold(1, 65535, "raw");
		run("Convert to Mask");
		run("Analyze Particles...", "  show=Overlay display");
		close("Results");
		close("Threshold");
		// Mask use
		open(mask);
		run("Set Scale...", "distance="+distance+" known="+known+" unit=µm");
		run("16-bit");
		run("Threshold...");
		setOption("BlackBackground", false);
		setThreshold(1, 65535, "raw");
		run("Convert to Mask");
		run("Analyze Particles...", "  show=Overlay display clear overlay");
		close("Results");
		close("Threshold");
		
		run("JACoP ", "imga="+imagename+"_IHC_mask.png imgb="+imagename+".png pearson");
		selectWindow("Log");
		localization_path = path+"results/colocalization/";
		if (!File.isDirectory(localization_path)) {
			File.makeDirectory(localization_path);
		}
		saveAs("Text", localization_path+imagename+".txt");
		close("Log");
		selectImage(imagename+"_IHC_mask.png");
		close();
		selectImage(imagename+".png");
		close();
	}
}

colocalization("64851", "3720", "27017.77");
colocalization("64852", "4320", "27071.77");
colocalization("64853", "3840", "27939.97");
colocalization("64854", "3840", "27939.97");
colocalization("64856", "3720", "27066.85");
colocalization("64857", "3960", "28813.10");
colocalization("16049-2-3101002", "3600", "26193.72");
colocalization("16049-3-3201002", "3720", "27066.85");
colocalization("16049-5-16110141B00201002", "3480", "25325.21");
colocalization("16049-6-3301002", "3720", "27066.85");
colocalization("16049-7-16110141B00301002", "3720", "27071.77");
colocalization("16049-8-3401002", "3480", "25320.60");
colocalization("16049-10-16110141B00401002", "3720", "27071.77");
colocalization("16045-1-101002", "4080", "29686.22");
colocalization("16045-2-16110141B00901002", "3840", "27945.06");
colocalization("16045-3-16110141B01001002", "3840", "27945.06");
colocalization("16045-4-201002", "3840", "27939.97");
colocalization("16045-5-16110141B01101002", "4440", "32311.47");
colocalization("16045-6-301002", "3600", "26193.72");
colocalization("16045-7-401003", "4320", "31432.47");
colocalization("16045-8-16110141B01201002", "3960", "28818.34");
colocalization("16045-9-16110141B01301002", "3312", "24091.65");
colocalization("16045-10-16110141B01401002", "3480", "25325.21");
colocalization("16045-11-16110141B01501002", "3480", "25325.21");
colocalization("16045-13-601002", "2720", "19780.02");
colocalization("16045-14-16110141B01601002", "3720", "27071.77");
colocalization("16045-7-701002", "4080", "29686.22");
colocalization("16045-17-801002", "3960", "28813.10");
colocalization("16045-18-901002", "4200", "30559.35");
colocalization("16045-19-1001002", "3960", "28813.10");
colocalization("16045-20-1101002", "3360", "24447.48");
colocalization("16045-21-1201002", "3240", "23574.35");
colocalization("16045-22-16110141B01701002", "2520", "18338.94");
colocalization("16045-23-16110141B01801002", "3240", "23578.64");
colocalization("16045-24-16110141B01901002", "3480", "25325.21");
colocalization("16045-25-1301002", "3480", "25320.60");
colocalization("16045-28-16110141B02001002", "3008", "21880.34");
colocalization("16045-29-16110141B02101002", "2880", "20958.79");
colocalization("16045-31-16110141B02201002", "3960", "28818.34");
colocalization("16045-32-1701002", "4080", "29686.22");
colocalization("16045-50-16110141B02401002", "3360", "24451.92");