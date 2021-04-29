import math
import os
from DICOMLib import DICOMUtils
import pydicom
from segmentation_dicom_python import segmentation_dicom

from __main__ import qt, slicer, vtk


#
# LungsImporter module
#

class LungsImporter:
  def __init__(self, parent):
    import string
    parent.title = "Lungs Importer"
    parent.categories = ["Experimental"]
    parent.contributors = ["NeuroLab"]
    parent.helpText = string.Template("""
This is experimental lungs DICOM and NII markup importer module
    """)
    parent.acknowledgementText = """
    Nothing now
    """
    self.parent = parent

#
# Widget
#


class LungsImporterWidget:

  def __init__(self, parent=None):
    self.parent = parent
    self.logic = None
    self.displayNode = None
    self.roiNode = None
    self.selectedDirectory=""
    self.series_list = []
    self.segments_list = []
    
  def setup(self):

    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )
    
    self.volumeNameSelectorLabel = qt.QLabel()
    self.volumeNameSelectorLabel.setText( "Select volume: " )
    self.volumeNameSelector = qt.QComboBox()
    self.volumeNameSelector.setEnabled(False)
    
    
    self.chooseDirectoryButton = qt.QPushButton("Choose directory")
    #self.analyzeDicomButton = qt.QPushButton("Analyze DICOM")
    self.importDicomAndNiiButton = qt.QPushButton("Import data")
    self.segmentBodyButton = qt.QPushButton("Segment body")
    
    self.lineEditPath = qt.QLineEdit()
    self.lineEditPath.setText("")
    self.labelPath = qt.QLabel()
    self.labelPath.setText("Path: ")
    
    self.segmentationProgressBar = qt.QProgressBar()
    self.segmentationProgressBar.setMinimum(0)
    self.segmentationProgressBar.setMaximum(100)
    self.segmentationProgressBar.setValue(0)
    self.segmentationProgressBar.setEnabled(False)
    
    layout.addRow(self.labelPath, self.lineEditPath) 
    layout.addRow(self.volumeNameSelectorLabel, self.volumeNameSelector) 
    layout.addRow(self.chooseDirectoryButton)
    layout.addRow(self.segmentBodyButton)
    layout.addRow(self.segmentationProgressBar)
    layout.addRow(self.importDicomAndNiiButton)
    
    self.chooseDirectoryButton.connect('clicked()', self.chooseDirectory)
    self.segmentBodyButton.connect('clicked()', self.segmentBody)
    self.importDicomAndNiiButton.connect('clicked()', self.importDicomAndNii)
    
  
  def segmentBody(self):
    thresh = -1800
    min_z = 10000
    ind = self.volumeNameSelector.currentIndex
    item = self.series_list[ind]
    volume_name = item["name"]
    
    if len(self.segments_list)<=0:
      niiFile = os.path.join(self.selectedDirectory, "segmentation.nii.gz")
      print("Nii will be: ", niiFile)
      segmentation_dicom.segmentation(self.selectedDirectory, niiFile, volume_name, self.segmentationProgressBar)
      self.segments_list.append(niiFile)
      print("Segmentation finished")
    else:
      self.segmentationProgressBar.setEnabled(True)
      self.segmentationProgressBar.setValue(100)
    
  
  def chooseDirectory(self):
    dialog = qt.QFileDialog(self.parent)
    filename = os.path.join(os.getcwd(),"lastpath.txt") 
    print(filename)
    try:
        f = open(filename, 'r')
        line = f.readline()
        f.close()
        print("Set dir", line)
        dialog.setDirectory(line.strip())
    except:
        print('Last path file not exists')
    dialog.setFileMode(qt.QFileDialog.Directory)
    if dialog.exec_():
        fnames = dialog.selectedFiles()
        if len(fnames)>=1:
            self.selectedDirectory = fnames[0]
    
    self.lineEditPath.setText(self.selectedDirectory)
    
    #Analyse for directory of CT and nii files
    #self.selectedDirectory = self.lineEditPath.text
       
    if not os.path.exists(self.selectedDirectory) or not os.path.isdir(self.selectedDirectory):
        print("Error, path not exists ", self.selectedDirectory)
        return
    
    f = open(filename, 'w')
    f.write("{}\n".format(self.selectedDirectory))
    f.close()
    
    files = os.listdir(self.selectedDirectory)
    
    self.dirCT = self.selectedDirectory
    #self.niiFile = os.path.join(self.selectedDirectory, "segmentation.nii.gz")
   
    segm_ext = "nii" 
    self.segments_list = []   
    for p in files:
        fp = os.path.join(self.selectedDirectory, p)
        
        sp = p.strip().split(".")
        if len(sp)>=2:
            if segm_ext in sp:
                self.segments_list.append(fp)
    
    if len(self.segments_list)>0:
      self.segmentationProgressBar.setEnabled(True)
      self.segmentationProgressBar.setValue(100)   
    else:
      self.segmentationProgressBar.setEnabled(True)
      self.segmentationProgressBar.setValue(0) 
            
    '''              
    if self.niiFile=="" or self.dirCT=="":
        print("Error, nii file or CT dir not detected in chosen dir files:", files)
        return
    else:
        print("CT dir:", self.dirCT)
        print("NII file:", self.niiFile)
    '''
     
    self.analyzeDICOMFile()
    
    
  def analyzeDICOMFile(self):
    dicomDataDir = self.dirCT  # input folder with DICOM files
    
    
    with DICOMUtils.TemporaryDICOMDatabase() as db:
      DICOMUtils.importDicom(dicomDataDir, db)
      patientUIDs = db.patients()
      
      uid = None
      
      for patientUID in patientUIDs:
        nm = db.nameForPatient(patientUID)
        print("Patient: ", nm)
        for study in db.studiesForPatient(patientUID):
            print("Study: ", study)
            series = db.seriesForStudy(study)
            print("Len = ", len(series))
            
            uid = series[0]
            
            for ser in series:
                files = db.filesForSeries(ser)
                dd = pydicom.dcmread(files[0])
                
                print(dd.SeriesDescription)
                
                seria_dict = {}
                seria_dict["name"] = dd.SeriesDescription
                seria_dict["uid"] = ser
                
                self.series_list.append(seria_dict)
    
    print(self.series_list)
                
    for item in self.series_list:
        self.volumeNameSelector.addItem(item["name"])
        
    self.volumeNameSelector.setCurrentIndex(0)
    self.volumeNameSelector.setEnabled(True)        


  def importDicomAndNii(self):

    dicomDataDir = self.dirCT  # input folder with DICOM files
        
        
    #Load segmentation  from .nii or .nii.gz files:
    for seg in self.segments_list:
      slicer.util.loadSegmentation(seg)
    
    ind = self.volumeNameSelector.currentIndex
    
    item = self.series_list[ind]
    
    with DICOMUtils.TemporaryDICOMDatabase() as db:
        DICOMUtils.importDicom(dicomDataDir, db)
        DICOMUtils.loadSeriesByUID([item['uid']])
