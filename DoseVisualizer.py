import math
import os
from DICOMLib import DICOMUtils
import pydicom
import numpy as np
import cv2
from scipy.io import savemat

from __main__ import qt, slicer, vtk

#
# LungsImporter module
#

class DoseVisualizer:
  def __init__(self, parent):
    import string
    parent.title = "Dose Visualizer"
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


class DoseVisualizerWidget:

  def __init__(self, parent=None):
    self.parent = parent
    self.logic = None
    self.displayNode = None
    self.roiNode = None
    self.selectedDirectory=""
    self.segmentationNode=None
    self.volumeNode=None
    
  def setup(self):

    print(qt)
    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )
    
    self.tableFrame = qt.QFrame()
    self.doseLayout = qt.QGridLayout()
    self.tableFrame.setLayout(self.doseLayout)
    self.parent.layout().addWidget(self.tableFrame)
    bins = [-1000, 0, 5, 10, 25, 50, 75, 100, 1000]
    self.color_widgets = []
    self.color_widgets.append(self.getColorWidgets("-1000", "0", 0.0, 0.0, 1.0))
    b = 0
    count = len(bins)-3
    step = 1.0/count
    g = 1.0
    r = 0.0
    for index in range(1, count+1):
        self.color_widgets.append(self.getColorWidgets(str(bins[index]), str(bins[index+1]), r, g, b))
        r+=step
        g-=step
        if r<0.0:
          r = 0.0
        if r>1.0:
          r = 1.0
        if g<0.0:
          g = 0.0
        if g>1.0:
          g = 1.0
          
    self.color_widgets.append(self.getColorWidgets("100", "1000", 1.0, 1.0, 1.0))
         
    self.mindlbl = qt.QLabel()
    self.mindlbl.setText("Min dose")
    self.maxdlbl = qt.QLabel()
    self.maxdlbl.setText("Max dose")
    self.clrlbl = qt.QLabel()
    self.clrlbl.setText("Area color")
    self.doseLayout.addWidget(self.mindlbl, 0, 0) 
    self.doseLayout.addWidget(self.maxdlbl, 0, 1) 
    self.doseLayout.addWidget(self.clrlbl, 0, 2) 
    row = 1
    for d in self.color_widgets:
      self.doseLayout.addWidget(d["min_lbl"], row, 0)
      self.doseLayout.addWidget(d["max_lbl"], row, 1)
      self.doseLayout.addWidget(d["line_ed"], row, 2)
      row+=1
    
    self.selectorLabel = qt.QLabel()
    self.selectorLabel.setText( "Reference volume: " )
    self.volumeSelector = slicer.qMRMLNodeComboBox()
    self.volumeSelector.nodeTypes = ( "vtkMRMLScalarVolumeNode", "" )
    self.volumeSelector.noneEnabled = False
    self.volumeSelector.selectNodeUponCreation = True
    self.volumeSelector.setMRMLScene( slicer.mrmlScene )
    self.volumeSelector.setToolTip( "Select reference volume for dose calc" )
       
    
    
    self.bodySegmentSelectorLabel = qt.QLabel()
    self.bodySegmentSelectorLabel.setText( "Select segmentation: " )
    self.bodySegmentSelector = slicer.qMRMLNodeComboBox()
    self.bodySegmentSelector.nodeTypes = ( "vtkMRMLSegmentationNode", "" )
    self.bodySegmentSelector.noneEnabled = False
    self.bodySegmentSelector.selectNodeUponCreation = True
    self.bodySegmentSelector.setMRMLScene( slicer.mrmlScene )
    
    self.doseSegmentSelectorLabel = qt.QLabel()
    self.doseSegmentSelectorLabel.setText( "Isodose segment: " )
    self.doseSegmentSelector = slicer.qMRMLNodeComboBox()
    self.doseSegmentSelector.nodeTypes = ( "vtkMRMLSegmentationNode", "" )
    self.doseSegmentSelector.noneEnabled = False
    self.doseSegmentSelector.selectNodeUponCreation = True
    self.doseSegmentSelector.setMRMLScene( slicer.mrmlScene )
    
    
    self.showRadiationMarkerCheckBox = qt.QCheckBox()
    self.showRadiationMarkerCheckBox.setText("Show test marker")
    
    self.visualizeSegmentButton = qt.QPushButton("Visualize segment")
    self.generateIsolinesButton = qt.QPushButton("Generate isolines")
    #self.visualizeIsodoseButton = qt.QPushButton("Visualize isodose segment")
    
    
    layout.addRow(self.selectorLabel, self.volumeSelector)  
    layout.addRow(self.bodySegmentSelectorLabel, self.bodySegmentSelector) 
    #layout.addRow(self.doseSegmentSelectorLabel, self.doseSegmentSelector) 
    layout.addRow(self.showRadiationMarkerCheckBox)
    layout.addRow(self.visualizeSegmentButton)
    layout.addRow(self.generateIsolinesButton)
    #layout.addRow(self.visualizeIsodoseButton)
    
    self.visualizeSegmentButton.connect('clicked()', self.visualizeSegment)
    self.generateIsolinesButton.connect('clicked()', self.generateIsolines)
    #self.visualizeIsodoseButton.connect('clicked()', self.visualizeIsodose)
    self.showRadiationMarkerCheckBox.connect('clicked(bool)', self.showRadiationMarkerCheckBoxChecked)
  
  def getColorWidgets(self, min_txt, max_txt, r,g,b):
    color = {}
    color['r']=r
    color['g']=g
    color['b']=b
    color["min_lbl"] = qt.QLabel()
    color["min_lbl"].setText(min_txt)
    color["max_lbl"] = qt.QLabel()
    color["max_lbl"].setText(max_txt)
    color["line_ed"] = qt.QLineEdit()
    color["line_ed"].setText("")
    color["line_ed"].setReadOnly(True)
    color["line_ed"].setStyleSheet("background: rgb({}, {}, {});".format(int(r*255), int(g*255), int(b*255))) 
    color["line_ed"].setFixedWidth(150)
    return color
  def visualizeIsodose(self):
    #self.segmentationNode = self.doseSegmentSelector.currentNode()
    self.segmentationNode.GetSegmentation().SetConversionParameter('Smoothing factor','1.0')
    segmentationDisplayNode = self.segmentationNode.GetDisplayNode()
    segmentationDisplayNode.SetAllSegmentsVisibility(True)
    visibleSegmentIds = vtk.vtkStringArray()
    self.segmentationNode.GetDisplayNode().GetVisibleSegmentIDs(visibleSegmentIds)
    '''
    for segmentIndex in range(visibleSegmentIds.GetNumberOfValues()):
        segmentID = visibleSegmentIds.GetValue(segmentIndex)
        if segmentIndex<1:
          segmentationDisplayNode.SetSegmentVisibility(segmentID, False)
    '''
    
    visibleSegmentIds = vtk.vtkStringArray()
    self.segmentationNode.GetDisplayNode().GetVisibleSegmentIDs(visibleSegmentIds)
    segments_ids = []
    b = 0
    count = visibleSegmentIds.GetNumberOfValues()
    step = 1.0/count
    g = 1.0
    r = 0.0
    for segmentIndex in range(visibleSegmentIds.GetNumberOfValues()):
        segmentID = visibleSegmentIds.GetValue(segmentIndex)
        segmentationDisplayNode.SetSegmentOverrideColor(segmentID, r, g, b)
        r+=step
        g-=step
        if r<0.0:
          r = 0.0
        if r>1.0:
          r = 1.0
        if g<0.0:
          g = 0.0
        if g>1.0:
          g = 1.0
        print(segmentID)
        
    segmentationDisplayNode.SetSegmentVisibility(visibleSegmentIds.GetValue(0), False)
    self.segmentationNode.CreateClosedSurfaceRepresentation()
    segmentationDisplayNode.SetOpacity3D(0.1)
    #self.createMarkerForRadiation()
  
  
  
  def visualizeSegment(self):
    segmentationNode = self.bodySegmentSelector.currentNode()
    segmentationNode.GetSegmentation().SetConversionParameter('Smoothing factor','1.0')
    segmentationDisplayNode = segmentationNode.GetDisplayNode()
    segmentationDisplayNode.SetSegmentVisibility('Segment_1', True)
    segmentationNode.CreateClosedSurfaceRepresentation()
    segmentationDisplayNode.SetOpacity3D(0.3)
    self.createMarkerForRadiation()
    
    
  def showRadiationMarkerCheckBoxChecked(self, checked):
    print(checked)
    if checked:
      self.radiationRoi.GetDisplayNode().SetVisibility(True)
    else:
      self.radiationRoi.GetDisplayNode().SetVisibility(False)
  
  def generateIsolines(self):
    center = [0,0,0]
    radius = [0,0,0]
    self.radiationRoi.GetControlPointWorldCoordinates(0, center)
    self.radiationRoi.GetControlPointWorldCoordinates(1, radius)
    '''data = {'center': center, 'radius': radius}
    print(data)'''
    
    if self.volumeNode!=None:
      slicer.mrmlScene.RemoveNode(self.volumeNode)
      self.volumeNode=None
    
    referenceVolume = self.volumeSelector.currentNode()
    
    imageSize = referenceVolume.GetImageData().GetDimensions()
    voxelType=vtk.VTK_UNSIGNED_CHAR
    imageOrigin = referenceVolume.GetOrigin()
    imageSpacing = referenceVolume.GetSpacing()
    imageDirections = [[0,0,0], [0,0,0], [0,0,0]]
    
    print("imdirs before", imageDirections)
    referenceVolume.GetIJKToRASDirections(imageDirections)
    print("imageDirections:", imageDirections)
    fillVoxelValue = 0
    
    self.radiationVolumeName = "RadiationVolume"
    
    # Create an empty image volume, filled with fillVoxelValue
    self.imageData = vtk.vtkImageData()
    self.imageData.SetDimensions(imageSize)
    self.imageData.AllocateScalars(voxelType, 1)
    self.thresholder = vtk.vtkImageThreshold()
    self.thresholder.SetInputData(self.imageData)
    self.thresholder.SetInValue(-1000)
    self.thresholder.SetOutValue(1000)
    self.thresholder.Update()
    # Create volume node
    self.volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", self.radiationVolumeName)
    self.volumeNode.SetOrigin(imageOrigin)
    self.volumeNode.SetSpacing(imageSpacing)
    self.volumeNode.SetIJKToRASDirections(imageDirections)
    self.volumeNode.SetAndObserveImageData(self.thresholder.GetOutput())
    self.volumeNode.CreateDefaultDisplayNodes()
    self.volumeNode.CreateDefaultStorageNode()
    
    '''
    mean_Ras = center
    print(mean_Ras)
    slicer.modules.markups.logic().JumpSlicesToLocation(mean_Ras[0], mean_Ras[1], mean_Ras[2], True)
    '''
    
    #Take the marker and try to render radiation volume
    corner1_Ras_global = [center[0]-radius[0], center[1]-radius[1], center[2]-radius[2]]
    corner2_Ras_global = [center[0]+radius[0], center[1]+radius[1], center[2]+radius[2]]
    
    # If volume node is transformed, apply that transform to get volume's RAS coordinates
    #transformRasToVolumeRas = vtk.vtkGeneralTransform()
    #slicer.vtkMRMLTransformNode.GetTransformBetweenNodes(None, self.volumeNode.GetParentTransformNode(), transformRasToVolumeRas)
    #corner1_Ras = transformRasToVolumeRas.TransformPoint(corner1_Ras_global[0:3])
    #corner2_Ras = transformRasToVolumeRas.TransformPoint(corner2_Ras_global[0:3])
        
    
    # Get voxel coordinates from physical coordinates
    volumeRasToIjk = vtk.vtkMatrix4x4()
    self.volumeNode.GetRASToIJKMatrix(volumeRasToIjk)
    
    corner1_Ijk = [0, 0, 0, 1]
    corner2_Ijk = [0, 0, 0, 1]
    
    volumeRasToIjk.MultiplyPoint(np.append(corner1_Ras_global,1.0), corner1_Ijk)
    volumeRasToIjk.MultiplyPoint(np.append(corner2_Ras_global,1.0), corner2_Ijk)
    
    point_corner1_Ijk = [ int(round(c)) for c in corner1_Ijk[0:3] ]
    point_corner2_Ijk = [ int(round(c)) for c in corner2_Ijk[0:3] ]
    
    volumeImageExtent = self.imageData.GetExtent()
    print(volumeImageExtent)
    
    point_corner1_Ijk_volume = [point_corner1_Ijk[0]-volumeImageExtent[0], point_corner1_Ijk[1]-volumeImageExtent[2], point_corner1_Ijk[2]-volumeImageExtent[4]]
    point_corner2_Ijk_volume = [point_corner2_Ijk[0]-volumeImageExtent[0], point_corner2_Ijk[1]-volumeImageExtent[2], point_corner2_Ijk[2]-volumeImageExtent[4]]
    
    point_corner1_Ijk_volume.reverse()
    point_corner2_Ijk_volume.reverse()
    
    d_x = abs(point_corner2_Ijk_volume[0] - point_corner1_Ijk_volume[0])
    d_y = abs(point_corner2_Ijk_volume[1] - point_corner1_Ijk_volume[1])
    d_z = abs(point_corner2_Ijk_volume[2] - point_corner1_Ijk_volume[2])
    
    #x_axis = np.linspace(-1, 1, d_x)[:, None, None]
    #y_axis = np.linspace(-1, 1, d_y)[None, :, None]
    #z_axis = np.linspace(-1, 1, d_z)[None, None, :]
    
    x_axis = np.abs(np.linspace(-256, 256, d_x))
    y_axis = np.abs(np.linspace(-256, 256, d_y))
    z_axis = np.abs(np.linspace(-256, 256, d_z))
    [X,Y,Z] = np.meshgrid(x_axis,y_axis,z_axis,indexing = 'ij');
    
    arr = np.sqrt(X ** 2 + Y ** 2 + Z**2)
    
    arr = arr/np.max(arr)
    
    #Inverse, large values inside
    
    arr = np.subtract(1.0, arr)
    arr = arr*100.0
    
    print("min", np.min(arr), "max", np.max(arr))
    
    #digitize an array (reverse quadratic)
    '''
    bins_temp = []
    for x in range(0, 11):
      bins_temp.append(x*x)
    bins_temp.reverse()
    bins = [100-x for x in bins_temp]
    '''
    bins = [-1000, 0, 5, 10, 25, 50, 75, 100, 1000]
    arr = np.digitize(arr, bins)
    
    
    #arr = arr.round()
    #print(arr[int(d_x/2), :, :])
    
    #cv2.imshow("arr", arr[int(d_x/2), :, :])
        
    voxelArray = slicer.util.arrayFromVolume(self.volumeNode)
    
    origin = self.volumeNode.GetOrigin()
    spacing = self.volumeNode.GetSpacing()
    
    print(origin, spacing)
    
    print(center, radius, point_corner1_Ijk_volume, point_corner2_Ijk_volume)
    #print(arr.shape)
    
    #print(voxelArray.shape)
    
    с1 = [min(point_corner1_Ijk_volume[0], point_corner2_Ijk_volume[0]), min(point_corner1_Ijk_volume[1], point_corner2_Ijk_volume[1]), min(point_corner1_Ijk_volume[2], point_corner2_Ijk_volume[2])]
    с2 = [max(point_corner1_Ijk_volume[0], point_corner2_Ijk_volume[0]), max(point_corner1_Ijk_volume[1], point_corner2_Ijk_volume[1]), max(point_corner1_Ijk_volume[2], point_corner2_Ijk_volume[2])]
    
    print(voxelArray[с1[0]:с2[0], с1[1]:с2[1], с1[2]:с2[2]].shape)
    voxelArray[с1[0]:с2[0], с1[1]:с2[1], с1[2]:с2[2]] = arr
    self.volumeNode.Modified()
    
    '''
    d = {"volume_array":arr}
    savemat("volume.mat", d)
    '''
    
    filename = os.path.join(os.getcwd(),"volume.nrrd") 
    
    myStorageNode = self.volumeNode.CreateDefaultStorageNode()
    myStorageNode.SetFileName(filename)
    myStorageNode.WriteData(self.volumeNode)
    
    if self.segmentationNode!=None:
      slicer.mrmlScene.RemoveNode(self.segmentationNode)
      self.segmentationNode=None
    
    self.segmentationNode = slicer.util.loadSegmentation(filename)
    self.segmentationNode.SetName("DoseSegmentation")
    self.visualizeIsodose()

    
    
  def createMarkerForRadiation(self):
    referenceVolume = self.volumeSelector.currentNode()
    if referenceVolume==None:
      return
      
    imageSize = referenceVolume.GetImageData().GetDimensions()
    imageOrigin = referenceVolume.GetOrigin()
    
    size = [int(float(x)/10.0) for x in imageSize]
    
    print(imageSize, imageOrigin)
    
    self.radiationRoi=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLAnnotationROINode")
    self.radiationRoi.SetName('Radiation ROI')
    self.radiationRoi.SetXYZ(imageOrigin[0], imageOrigin[1], imageOrigin[2])
    self.radiationRoi.SetRadiusXYZ(size)
    #self.radiationRoi.GetDisplayNode().SetVisibility(False)
    
    
