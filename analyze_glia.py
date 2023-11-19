# Microglia Quantification v9.2
#
# This ImageJ script measures microglia in images. 
#
# Usage:
# See Installation below. To run this script, drag-and-drop the script file to Fiji. 
# Then, run the script using Ctrl+R or from the menu.
# It will open a window and ask the user to pick a folder (of nd2/tif images). 
# 
# Note:
# Each image must have 2 channels: channel 1 DAPI, channel 2 microglia staining.
# For every image the script performs theses steps:
#  - Maximal Intensity Projection from Z dimension.
#  - Find nuclei using channel 1.
#  - Find microglia using channel 2 and verify that there is a nucleus.
#  - Measures microglia area, perimeter, RI, cable-length. 
#    - Area [um^2] of cell (after Maximum Intensity Projection)
#    - Perimeter [um] - the 2D length of the memebrane of the cell
#    - RI = (perimeter/area)/[2(π/area)^(1/2)]
#    - Cable-length is the sum of lengths of a cell's skeleton ramifocations.
#  - Export to excel file in the folder.
#
# Installation:
#  - Copy membrane.ilp to "C:\Users\myusername\microglia\" , where "myusername" is your user name.
#  - Install the Fiji plugins listed below.
#
# Fiji plugins requirements:
#  - ResultsToExcel (for Read and Write Excel)
#  - Neuroanatomy (for SNT)
#  - ROIs to Masks
#  - CLIJ and CLIJ2
#  - ilastik (and "Configure ilastik executable location")
#  - Masks from ROIs
#
# Changes in v2: 
#  - Fixed cable_length in excel.
#  - Fixed open files with spaces.
# Changes in v3: 
#  - Increased background subtraction and also to DAPI.
#  - Create folder summary in Excel file
# Changes in v4: 
#  - Fixed: get_cable_length() no longer raises NullPointerException, and returns 0
#  - Fixed: File not found .csv when creating folder summary.
#  - Cleaned imports.
# Changes in v6: 
#  - support .tif
# Changes in v8: 
#  - Added Ilastik membrane mask to segment ramifications
#  - Added MICROGLIA_MAX_SIZE

import os
from math import pi as PI, sqrt
from ij import IJ, ImagePlus, Prefs
from ij.process import ImageConverter
from ij import WindowManager
from ij.plugin.frame import RoiManager
from ij.measure import Measurements, ResultsTable
from sc.fiji.snt.analysis import SkeletonConverter, TreeAnalyzer, SNTTable
from sc.fiji.snt import Tree
from ij.gui import GenericDialog, WaitForUserDialog


DEBUG = False
# DEBUG = True
CONNECT_MICROGLIA = 0.8  # Length [microns] of gaps to connect when segmenting microglia
NUCLEUS_INTERSECTION_WITH_CELL = 15  # Minimum size of nucleus in cell [microns]
MICROGLIA_MIN_SIZE = 20
MICROGLIA_MIN_SIZE_3D = 10
MICROGLIA_MAX_SIZE = 500
MEMBRANE_MODEL_FILE =  "~\\microglia\\membrane.ilp"  # ~ exapnds to the user's dicrectory, e.g.: C:\Users\myusername


def get_cable_length(roi, ip, imp):
    ip.setRoi(roi);
    cal = imp.getCalibration()
    cropped = ImagePlus("cropped image", ip.crop())
    cropped.setCalibration(cal)
    # cropped.show()
    # IJ.showMessage("1111")
    skelConv = SkeletonConverter(cropped, True)
    # Debug: IJ.log(str(cropped.getRawStatistics().max))
    cable_length = 0.0
    if cropped.getRawStatistics().max > 0:
        #    IJ.log(str(skelConv))
        #    IJ.log(str(dir(skelConv)))
        trees = skelConv.getTrees()
        #  print "trees", trees
        tree = trees[0]
        if len(trees) > 1:
            # cropped.show()
            # print "tree lengths", [len(x.list())  for x in trees]
            max_tree_len = max(len(x.list())  for x in trees)
            # Find largest tree in ROI
            tree = next(x for x in trees if max_tree_len == len(x.list()))
        # print "longest tree", tree
        analyzer = TreeAnalyzer(tree)
        cable_length = analyzer.getMetric("Cable length")
    cropped.close()
    return cable_length


def remove_inverted_lut():
    imp = IJ.getImage()
    ip = imp.getProcessor()
    print "isInvertedLut", imp.isInvertedLut()
    if imp.isInvertedLut():
        print "inverting"
        IJ.run(imp, "Invert LUT", "")


def main():
    # Reset
    IJ.run("Close All");
    IJ.run("Set Measurements...", "area area_fraction perimeter display redirect=None decimal=3");
    RM = RoiManager()        # we create an instance of the RoiManager class
    rm = RM.getRoiManager()
    rm.reset();
    
    # Choose folder
    input_dir = IJ.getDirectory("Select a folder with images");
    print "input_dir: ", input_dir

    # Loop images
    file_list = os.listdir(input_dir)
    print file_list
    print "Found ", len(file_list), "files"
    analysis_dir = os.path.join(input_dir, "analysis_v9")
    if not os.path.exists(analysis_dir):
        os.mkdir(analysis_dir)  
    
    for file_name in file_list:
        file_path = os.path.join(input_dir, file_name); 
        print("loop " + file_path);
        if not (file_name.endswith(".nd2") or file_name.endswith(".tif")):
          continue
        
        rm.reset();  # Reset ROIs in Manager
        # Open image
        print "processing image: ", file_path
        IJ.run("Bio-Formats Importer", 
               "open='" + file_path + 
               "' autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT")
        if WindowManager.getImageCount() == 0:
            print "Could not open image " + file_path
            IJ.showMessage("Could not open image " + file_path)
            exit()
        IJ.run("Z Project...", "projection=[Max Intensity]");
        
        # MIP
        imp = IJ.getImage()
        imp.setTitle("MIP");
        
        # Nuceli
        IJ.run("Duplicate...", "duplicate channels=1");
        imp = IJ.getImage()
        imp.setTitle("nuclei")
        IJ.run("Subtract Background...", "rolling=30")
        IJ.setAutoThreshold(imp, "Otsu dark")
        IJ.run(imp, "Convert to Mask", "")
        remove_inverted_lut()

        # Microglia segmentation
        IJ.selectWindow("MIP")
        IJ.run("Duplicate...", "duplicate channels=2-2")
        imp = IJ.getImage()
        imp.setTitle("antibody_raw")
        calibration = imp.getCalibration()
        print("calibration", calibration)
        IJ.run("Enhance Contrast", "saturated=0.35")
                
        membrane_model_file = os.path.expanduser(MEMBRANE_MODEL_FILE)
        assert os.path.exists(membrane_model_file), membrane_model_file
        if not DEBUG:
            IJ.run("Run Pixel Classification Prediction", 
                        "projectfilename=[" + membrane_model_file + "] " + 
                        " inputimage=[" + file_path + "] pixelclassificationtype=Segmentation")
        else:
            segmented_file = os.path.join(analysis_dir, file_name + "_seg.tif")
            if os.path.exists(segmented_file):
                IJ.open(segmented_file)
            else:
                IJ.run("Run Pixel Classification Prediction", 
                        "projectfilename=" + membrane_model_file + " " + 
                        " inputimage=[" + file_path + "] pixelclassificationtype=Segmentation")
                IJ.save(segmented_file)
        IJ.getImage().setTitle("segmented_3D")  # 1 = backgound, 2 = foreground cell
        imp = IJ.getImage()
        IJ.setRawThreshold(imp, 2, 255);
        IJ.run(imp, "Convert to Mask", "method=Default background=Dark ");
        IJ.run("Duplicate...", "duplicate")
        IJ.getImage().setTitle("objects_3D")
        IJ.run(IJ.getImage(), "Invert", "stack")
        IJ.run("3D OC Options", " redirect_to=none")
        IJ.run("3D Objects Counter on GPU (CLIJx, Experimental)", "cl_device=[Quadro M4000] threshold=254 slice=4 " + 
                " min.=" + str(MICROGLIA_MIN_SIZE_3D) + " max.=9999999 objects");
        IJ.run("Enhance Contrast", "saturated=0.35");
        IJ.run("Z Project...", "projection=[Max Intensity]");
        IJ.run(IJ.getImage(), "Fire", "");
        ImageConverter.setDoScaling(False);
        IJ.run("16-bit")
        IJ.setRawThreshold(IJ.getImage(), 1, 65536)
        IJ.run(IJ.getImage(), "Convert to Mask", "")
        imp_objects_2d = IJ.getImage()
        imp_objects_2d.setCalibration(calibration)
        print("calibration", calibration)
        IJ.run("Analyze Particles...", "size=" + 
               str(MICROGLIA_MIN_SIZE) + "-" + str(MICROGLIA_MAX_SIZE) + " show=Masks exclude add slice")
        IJ.getImage().setTitle("microglia_antibody_mask")
        remove_inverted_lut()

        # Remove microglia ROIs that do not have a significant nucleus
        IJ.selectWindow("nuclei");
        rm.runCommand("Show All with labels")
        imp = IJ.getImage()
        to_be_deleted = []
        print rm.getCount()
        if rm.getCount() > 0:
            for i in range(rm.getCount()):
                roi = rm.getRoi(i)
                imp.setRoi(roi)
                stats = imp.getAllStatistics();  # getStatistics returns are of ROI and not foreground area
                if DEBUG:
                    # print dir(stats)
                    print "areaFraction", stats.areaFraction
                sum_nuclei_gray = stats.area * stats.areaFraction / 100  # area of nucleus in cell calibrated
                if DEBUG:
                    print "sum_nuclei_gray", sum_nuclei_gray
                if sum_nuclei_gray < NUCLEUS_INTERSECTION_WITH_CELL:
                    to_be_deleted.append(i)
            if DEBUG:
                print "to_be_deleted", to_be_deleted
            rm.setSelectedIndexes(to_be_deleted)
            rm.runCommand("Delete");
            
        # User interaction
        myWait = WaitForUserDialog("Microglia Segmentation", "You can now edit the automatically detected segments and click OK.")
        myWait.show()

        # Measure Microglia ROIs
        IJ.selectWindow("microglia_antibody_mask")
        if rm.getCount() > 0:
            # Convert ROI manager to Mask
            IJ.run("Select All")
            IJ.run("Clear", "slice")
            # Create new window with ROIs as mask:
            IJ.run("Binary (0-255) mask(s) from Roi(s)", "show_mask(s) save_in=[] suffix=[] rm=[RoiManager[visible=true]]")
        else:
            IJ.run("Duplicate...", "duplicate ignore ")
        imp = IJ.getImage()  # without this setTitle doesn't work
        imp.setTitle("microglia_mask")
        if rm.getCount() > 0:
            # Measure microglia
            IJ.selectWindow("microglia_mask")
            IJ.run("Clear Results")
            rm.runCommand("Measure");
            # Add RI measure
            results = ResultsTable.getResultsTable()
            print "results.size()", results.size()
            imp = IJ.getImage()
            ip = imp.getProcessor()
            for row in range(results.size()):
                perimeter = results.getValue("Perim.", row)
                area = results.getValue("Area", row)  # ROI area
                ri = perimeter / area / (2 * sqrt((PI / area)))
                assert ri >= 1, ri  # Assert. Should never happen.
                results.setValue("RI", row, ri)
                assert results.size() == rm.getCount()  # to make sure that ROIs match results
                roi = rm.getRoi(row)
                IJ.log(str(row)+ " " + str(roi));
                cable_length = get_cable_length(roi, ip, imp)
                results.setValue("Cable_Length", row, cable_length)
            results.show("Results")
            results.deleteColumn("%Area")
            
            # Save CSV
            csv_file = os.path.join(analysis_dir, file_name + ".csv")
            print "Saving csv_file " + csv_file
            IJ.saveAs("Results", csv_file);
            
            # Save excel
            excel_file = os.path.join(analysis_dir, file_name + ".xlsx")
            IJ.run("Summarize");
            print "Saving image results " + excel_file
            IJ.run("Read and Write Excel", "file=[" + excel_file + "] dataset_label=[]")

        # Save image with segmentation.
        IJ.selectWindow("MIP");
        IJ.run("Make Composite");
        IJ.getImage().setC(1)  # set channel
        IJ.run("Enhance Contrast", "saturated=0.35")
        IJ.getImage().setC(2)
        IJ.run("Enhance Contrast", "saturated=0.35")
        
        IJ.run("Show Overlay");
        mip_image_file = analysis_dir + "/" + file_name + ".tif"
        print "Saving " + mip_image_file
        IJ.save(mip_image_file);
        
        IJ.selectWindow("microglia_mask")
        IJ.run("Merge Channels...", "c1=microglia_mask c3=nuclei create keep")
        IJ.save(analysis_dir + "/" + file_name + "_masks.png")
        # return
        IJ.run("Close All")  # close windows after processing each image
        # end processing one image
    
    # Create folder Summary
    IJ.run("Clear Results")
    results = ResultsTable.getResultsTable()
    for file_name in file_list:
        if not (file_name.endswith(".nd2") or file_name.endswith(".tif")):
          continue
        file_path = os.path.join(input_dir, file_name) 
        csv_file = os.path.join(analysis_dir, file_name + ".csv")
        print csv_file
        if not os.path.isfile(csv_file):
          continue
        IJ.open(csv_file)
        csv_table = ResultsTable.getActiveTable()
        print "heads", csv_table.headings
        for row in range(csv_table.size()):
            results.incrementCounter()
            results.addLabel(str(csv_table.getLabel(row)))
            for col in range(1, len(csv_table.headings)):
                val = csv_table.getValue(csv_table.headings[col], row)
                results.addValue(csv_table.headings[col], val)
        IJ.selectWindow(csv_table.getTitle()); 
        IJ.run("Close");
    results.show("Results")
    IJ.run("Summarize");
    
    # parent_name = os.path.normpath(input_dir).split(os.path.sep)[-1]
    # print "parent_name", parent_name
    # folder_excel_file = os.path.join(analysis_dir, "summary_"  + 
    #     parent_name + ".csv")
    # print "Saving folder results " + folder_excel_file
    # if not os.path.exists(analysis_dir):
    #         os.mkdir(analysis_dir)  
    # IJ.saveAs("Results", folder_excel_file);    
    # Theres is a bug in "Read and Write Excel" that doesn't export labels
    # IJ.run("Read and Write Excel", "file=[" + folder_excel_file + "] dataset_label=[]")
        
    
    print "Finished"
    IJ.showMessage("Finshed processing folder.")
    # end loop on images

main()
