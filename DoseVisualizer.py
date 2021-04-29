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
    
  def setup(self):

    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )
    
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
    
    
    self.showRadiationMarkerCheckBox = qt.QCheckBox()
    self.showRadiationMarkerCheckBox.setText("Show test marker")
    
    self.visualizeSegmentButton = qt.QPushButton("Visualize segment")
    self.generateIsolinesButton = qt.QPushButton("Generate isolines")
    
    
    layout.addRow(self.selectorLabel, self.volumeSelector)  
    layout.addRow(self.bodySegmentSelectorLabel, self.bodySegmentSelector) 
    layout.addRow(self.showRadiationMarkerCheckBox)
    layout.addRow(self.visualizeSegmentButton)
    layout.addRow(self.generateIsolinesButton)
    
    self.visualizeSegmentButton.connect('clicked()', self.visualizeSegment)
    self.generateIsolinesButton.connect('clicked()', self.generateIsolines)
    self.showRadiationMarkerCheckBox.connect('clicked(bool)', self.showRadiationMarkerCheckBoxChecked)
  
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
    
    #arr = np.subtract(1.0, arr)
    #arr = arr*256
    #arr = arr.round()
    #print(arr[int(d_x/2), :, :])
    
    cv2.imshow("arr", arr[int(d_x/2), :, :])
        
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
    
    d = {"volume_array":arr}
    savemat("volume.mat", d)

    
    
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
    
    
