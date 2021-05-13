import math
import os
from DICOMLib import DICOMUtils
import pydicom
from segmentation_dicom_python import segmentation_dicom

from __main__ import qt, slicer, vtk


#
# LungsImporter module
#

class Lungs3DCommon:
  def __init__(self, parent):
    import string
    parent.title = "Lungs 3D Common"
    parent.categories = ["Experimental"]
    parent.contributors = ["NeuroLab"]
    parent.helpText = string.Template("""
This is common united module for CT scan import and visualisation
    """)
    parent.acknowledgementText = """
    Nothing now
    """
    self.parent = parent

#
# Widget
#


class Lungs3DCommonWidget:

  def __init__(self, parent=None):
    self.parent = parent
    self.logic = None
    self.selectedDirectory=""
    self.series_list = []
    self.segments_list = []
    self.displayNode = None
    self.roiNode = None
    self.maskedVolumeNode = None
    
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
    
    self.segmentLabel = qt.QLabel()
    self.segmentLabel.setText( "Select segment: " )
    self.segmentSelector = slicer.qMRMLNodeComboBox()
    self.segmentSelector.nodeTypes = ( "vtkMRMLSegmentationNode", "" )
    self.segmentSelector.noneEnabled = False
    self.segmentSelector.selectNodeUponCreation = True
    self.segmentSelector.setMRMLScene( slicer.mrmlScene )
    #self.segmentSelector.setToolTip( "Pick the markup list to be filled" )
    layout.addRow(self.segmentLabel, self.segmentSelector)
    
    self.selectorLabel = qt.QLabel()
    self.selectorLabel.setText( "Select volume: " )
    self.volumeSelector = slicer.qMRMLNodeComboBox()
    self.volumeSelector.nodeTypes = ( "vtkMRMLScalarVolumeNode", "" )
    self.volumeSelector.noneEnabled = False
    self.volumeSelector.selectNodeUponCreation = True
    self.volumeSelector.setMRMLScene( slicer.mrmlScene )
    self.volumeSelector.setToolTip( "Pick the markup list to be filled" )
    layout.addRow(self.selectorLabel, self.volumeSelector) 
    
    self.ctTypeSelectorLabel = qt.QLabel()
    self.ctTypeSelectorLabel.setText( "Select CT type: " )
    self.ctTypeSelector = qt.QComboBox()
    self.ctTypeSelector.addItem('CT-Muscle')
    self.ctTypeSelector.addItem('CT-Lungs')
    self.ctTypeSelector.addItem('CT-Cardiac')
    self.ctTypeSelector.addItem('CT-Cardiac2')
    self.ctTypeSelector.addItem('CT-Coronary-Arteries-2')
    self.ctTypeSelector.setCurrentIndex(0)
    layout.addRow(self.ctTypeSelectorLabel, self.ctTypeSelector) 
    
    self.renderVolumeButton = qt.QPushButton("Render masked volume")
    self.calcTubePositionButton = qt.QPushButton("Calc tube position")
    layout.addRow(self.renderVolumeButton)
    layout.addRow(self.calcTubePositionButton) 
    
    self.chooseDirectoryButton.connect('clicked()', self.chooseDirectory)
    self.segmentBodyButton.connect('clicked()', self.segmentBody)
    self.importDicomAndNiiButton.connect('clicked()', self.importDicomAndNii)
    
    self.renderVolumeButton.connect('clicked()', self.renderVolume)
    self.calcTubePositionButton.connect('clicked()', self.calcTubePosition)

############################################################  
  
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
      self.segmentationProgressBar.setValue(100)
      print("Segmentation finished")
    else:
      self.segmentationProgressBar.setEnabled(True)
      self.segmentationProgressBar.setValue(100)

############################################################   
  
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
  
############################################################  
    
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

############################################################

  def importDicomAndNii(self):

    dicomDataDir = self.dirCT  # input folder with DICOM files
        
        
    #Load segmentation  from .nii or .nii.gz files:
    for seg in self.segments_list:
      seg_node = slicer.util.loadSegmentation(seg)
      self.segmentSelector.setCurrentNode(seg_node)
    
    ind = self.volumeNameSelector.currentIndex
    
    item = self.series_list[ind]
    
    with DICOMUtils.TemporaryDICOMDatabase() as db:
        DICOMUtils.importDicom(dicomDataDir, db)
        volume_nodes = DICOMUtils.loadSeriesByUID([item['uid']])
        #print("volume_node_name", volume_node_name)
        volume_node = slicer.mrmlScene.GetNodeByID(volume_nodes[0])
        #print("volume_node", volume_node)
        self.volumeSelector.setCurrentNode(volume_node)
        
############################################################
   
  def renderVolume(self):
    segmentationNode = self.segmentSelector.currentNode()  
    self.volumeNode = self.volumeSelector.currentNode()
    
    #Show segment
    segmentationDisplayNode = segmentationNode.GetDisplayNode()
    segmentationDisplayNode.SetSegmentVisibility('Segment_1', True)
    segmentationDisplayNode.SetSegmentVisibility('Segment_2', False)

    # Write segmentation to labelmap volume node with a geometry that matches the volume node
    labelmapVolumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode')
    slicer.modules.segmentations.logic().ExportVisibleSegmentsToLabelmapNode(segmentationNode, labelmapVolumeNode, self.volumeNode)

    # Masking
    import numpy as np
    voxels = slicer.util.arrayFromVolume(self.volumeNode)
    mask = slicer.util.arrayFromVolume(labelmapVolumeNode)
    maskedVoxels = np.copy(voxels)  # we don't want to modify the original volume
    maskedVoxels[mask==0] = -1000

    # Write masked volume to volume node and show it
    if self.maskedVolumeNode!=None:
        slicer.mrmlScene.RemoveNode(self.maskedVolumeNode)
        self.maskedVolumeNode=None
        
    self.maskedVolumeNode = slicer.modules.volumes.logic().CloneVolume(self.volumeNode, "Masked")
    slicer.util.updateVolumeFromArray(self.maskedVolumeNode, maskedVoxels)
    slicer.util.setSliceViewerLayers(self.maskedVolumeNode)
    
    
    
    logic = slicer.modules.volumerendering.logic()
    #self.volumeNode = self.volumeSelector.currentNode()
    self.volumePreset = logic.GetPresetByName(self.ctTypeSelector.currentText)
    
    if self.roiNode == None:
        self.roiNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLAnnotationROINode")
        self.roiNode.SetName("Masked render ROI")
        
    if self.displayNode!=None:
        slicer.mrmlScene.RemoveNode(self.displayNode)
        self.displayNode=None
    
    if self.maskedVolumeNode!=None:
        self.displayNode = logic.CreateVolumeRenderingDisplayNode()
        self.displayNode.SetAndObserveROINodeID(self.roiNode.GetID())
        self.displayNode.SetCroppingEnabled(1)
        self.displayNode.UnRegister(logic)
        #print(self.displayNode)        
        slicer.mrmlScene.AddNode(self.displayNode)
        
        self.maskedVolumeNode.AddAndObserveDisplayNodeID(self.displayNode.GetID())
        logic.UpdateDisplayNodeFromVolumeNode(self.displayNode, self.maskedVolumeNode)
        self.propertyNode = self.displayNode.GetVolumePropertyNode()
        self.propertyNode.Copy(self.volumePreset)
        self.displayNode.SetVisibility(1)
        
############################################################

  def calcTubePosition(self):
    #Switch visibility
    segmentationNode = self.segmentSelector.currentNode() 
    segmentationDisplayNode = segmentationNode.GetDisplayNode()
    segmentationDisplayNode.SetSegmentVisibility('Segment_1', False)
    segmentationDisplayNode.SetSegmentVisibility('Segment_2', True)
    #Calc ROI
    self.calcSegmentPosition()
    
    #Render tube
    segmentationNode.GetSegmentation().SetConversionParameter('Smoothing factor','0.0')
    segmentationNode.CreateClosedSurfaceRepresentation()
  
############################################################
  
  def calcSegmentPosition(self):
    print("Fnc")
    segmentationNode = self.segmentSelector.currentNode() 

    # Compute bounding boxes
    import SegmentStatistics
    segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
    segStatLogic.getParameterNode().SetParameter("Segmentation", segmentationNode.GetID())
    segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.obb_origin_ras.enabled",str(True))
    segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.obb_diameter_mm.enabled",str(True))
    segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.obb_direction_ras_x.enabled",str(True))
    segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.obb_direction_ras_y.enabled",str(True))
    segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.obb_direction_ras_z.enabled",str(True))
    segStatLogic.computeStatistics()
    stats = segStatLogic.getStatistics()

    print("Statc calcd")

    # Draw ROI for each oriented bounding box
    import numpy as np
    for segmentId in stats['SegmentIDs']:
        print(segmentId)
        # Get bounding box
        obb_origin_ras = np.array(stats[segmentId,"LabelmapSegmentStatisticsPlugin.obb_origin_ras"])
        obb_diameter_mm = np.array(stats[segmentId,"LabelmapSegmentStatisticsPlugin.obb_diameter_mm"])
        obb_direction_ras_x = np.array(stats[segmentId,"LabelmapSegmentStatisticsPlugin.obb_direction_ras_x"])
        obb_direction_ras_y = np.array(stats[segmentId,"LabelmapSegmentStatisticsPlugin.obb_direction_ras_y"])
        obb_direction_ras_z = np.array(stats[segmentId,"LabelmapSegmentStatisticsPlugin.obb_direction_ras_z"])
        # get voxel size
        spacing = self.volumeNode.GetSpacing()
        print(spacing)    
        # fix for dicretization errors
        # adding ~1 voxel size offset to the end position
        obb_diameter_mm[2] = obb_diameter_mm[2] - spacing[1];

        # Create ROI
        segment = segmentationNode.GetSegmentation().GetSegment(segmentId)
        roi=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLAnnotationROINode")
        roi.SetName(segment.GetName()+' bounding box')
        roi.SetXYZ(0.0, 0.0, 0.0)
        roi.SetRadiusXYZ(*(0.5*obb_diameter_mm))
        
        # Position and orient ROI using a transform
        obb_center_ras = obb_origin_ras+0.5*(obb_diameter_mm[0] * obb_direction_ras_x + obb_diameter_mm[1] * obb_direction_ras_y + obb_diameter_mm[2] * obb_direction_ras_z)
        obb_end1center_ras = obb_center_ras - 0.5* obb_diameter_mm[2] * obb_direction_ras_z;
        obb_end2center_ras = obb_center_ras + 0.5* obb_diameter_mm[2] * obb_direction_ras_z;
        print('end1 end2 orientation')
        print(obb_end1center_ras)
        print(obb_end2center_ras)
        print(np.column_stack((obb_direction_ras_x, obb_direction_ras_y, obb_direction_ras_z)))

        boundingBoxToRasTransform = np.row_stack((np.column_stack((obb_direction_ras_x, obb_direction_ras_y, obb_direction_ras_z, obb_center_ras)), (0, 0, 0, 1)))
        boundingBoxToRasTransformMatrix = slicer.util.vtkMatrixFromArray(boundingBoxToRasTransform)
        transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
        transformNode.SetAndObserveMatrixTransformToParent(boundingBoxToRasTransformMatrix)
        roi.SetAndObserveTransformNodeID(transformNode.GetID())

        # Create new ROI
        roi=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode")
        roi.SetName(segment.GetName()+' bounding box')
        roi.SetXYZ(0.0, 0.0, 0.0)
        roi.SetRadiusXYZ(*(0.5*obb_diameter_mm))

        # construct cylinder coord system
        # here Z is pointing outward from the cylinder
        cyl_direction_ras_z = -1*np.sign(obb_direction_ras_z[2])*obb_direction_ras_z;
        cyl_direction_ras_z *= 1/np.linalg.norm(cyl_direction_ras_z)
        cyl_direction_ras_y = np.cross(cyl_direction_ras_z, np.append(cyl_direction_ras_z[0:2],0))
        cyl_direction_ras_y *= 1/np.linalg.norm(cyl_direction_ras_y)
        cyl_direction_ras_x = np.cross(cyl_direction_ras_y,cyl_direction_ras_z)
        cyl_direction_ras_x *= 1/np.linalg.norm(cyl_direction_ras_x)

        print('orientation2')
        print(np.column_stack((cyl_direction_ras_x, cyl_direction_ras_y, cyl_direction_ras_z)))
        
        obb_end1center_ras = obb_center_ras - 0.5* obb_diameter_mm[2] * cyl_direction_ras_z;
        obb_end2center_ras = obb_center_ras + 0.5* obb_diameter_mm[2] * cyl_direction_ras_z;
        boundingBoxToRasTransform = np.row_stack((np.column_stack((cyl_direction_ras_x, cyl_direction_ras_y, cyl_direction_ras_z, obb_center_ras)), (0, 0, 0, 1)))
        boundingBoxToRasTransformMatrix = slicer.util.vtkMatrixFromArray(boundingBoxToRasTransform)
        transformNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
        transformNode.SetAndObserveMatrixTransformToParent(boundingBoxToRasTransformMatrix)
        roi.SetAndObserveTransformNodeID(transformNode.GetID())
