import math
from __main__ import qt, slicer, vtk

#
# VolumetricExperiments module
#

class VolumetricExperiments:
  def __init__(self, parent):
    import string
    parent.title = "Volumetric Experiments"
    parent.categories = ["Experimental"]
    parent.contributors = ["NeuroLab"]
    parent.helpText = string.Template("""
This is experimental volumetric modelling control module
    """)
    parent.acknowledgementText = """
    Nothing now
    """
    self.parent = parent

#
# Widget
#


class VolumetricExperimentsWidget:

  def __init__(self, parent=None):
    self.parent = parent
    self.logic = None
    self.displayNode = None
    self.roiNode = None
    self.maskedVolumeNode = None
    
  def setup(self):

    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )
    
    # Markup selector
    #self.createVolumeButton = qt.QPushButton("Create test volume")
    #self.renderVolumeButton = qt.QPushButton("Render test volume")
    
    self.applySegmentButton = qt.QPushButton("Apply segment")
    #layout.addRow(self.createVolumeButton, self.renderVolumeButton)
    
    
    self.segmentLabel = qt.QLabel()
    self.segmentLabel.setText( "Select segment: " )
    self.segmentSelector = slicer.qMRMLNodeComboBox()
    self.segmentSelector.nodeTypes = ( "vtkMRMLSegmentationNode", "" )
    self.segmentSelector.noneEnabled = False
    self.segmentSelector.selectNodeUponCreation = True
    self.segmentSelector.setMRMLScene( slicer.mrmlScene )
    #self.segmentSelector.setToolTip( "Pick the markup list to be filled" )
    layout.addRow(self.segmentLabel, self.segmentSelector)
    
    '''
    self.tubeLabel = qt.QLabel()
    self.tubeLabel.setText( "Select tube segment: " )
    self.tubeSelector = slicer.qMRMLNodeComboBox()
    self.tubeSelector.nodeTypes = ( "vtkMRMLSegmentationNode", "" )
    self.tubeSelector.noneEnabled = False
    self.tubeSelector.selectNodeUponCreation = True
    self.tubeSelector.setMRMLScene( slicer.mrmlScene )
    #self.segmentSelector.setToolTip( "Pick the markup list to be filled" )
    layout.addRow(self.tubeLabel, self.tubeSelector)  
    '''
    
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
    
    '''
    self.lineEditMin = qt.QLineEdit()
    self.lineEditMin.setText(0)
    self.labelMin = qt.QLabel()
    self.labelMin.setText("Min thr: ")
    self.minButton = qt.QPushButton("Set min/max values")
    layout.addRow(self.labelMin, self.lineEditMin) 
    
    self.lineEditMax = qt.QLineEdit()
    self.lineEditMax.setText(0)
    self.labelMax = qt.QLabel()
    self.labelMax.setText("Max thr: ")
    #self.maxButton = qt.QPushButton("Set max")
    layout.addRow(self.labelMax, self.lineEditMax) 
    layout.addRow(self.minButton) 
    '''
    
    #self.calcTubeDataButton = qt.QPushButton("Calc tube data")
    
    #self.testSegmentROIParamsButton = qt.QPushButton("Test segments parameters")
    
    self.renderVolumeButton = qt.QPushButton("Render masked volume")
    self.calcTubePositionButton = qt.QPushButton("Calc tube position")
    
    #layout.addRow(self.applySegmentButton)
    #layout.addRow(self.renderVolumeButton)
    #layout.addRow(self.calcTubeDataButton)
    #layout.addRow(self.testSegmentROIParamsButton)
    #layout.addRow(self.applyRoiButton)

    layout.addRow(self.renderVolumeButton)
    layout.addRow(self.calcTubePositionButton)    
        
    #connections
    #self.createVolumeButton.connect('clicked()', self.createVolume)
    #self.renderVolumeButton.connect('clicked()', self.renderVolume)
    #self.applySegmentButton.connect('clicked()', self.applySegment)
    #self.calcTubeDataButton.connect('clicked()', self.calcTubeData)
    #self.testSegmentROIParamsButton.connect('clicked()', self.testSegmentROIParams)
    #self.applyRoiButton.connect('clicked()', self.applyROI)
    
    self.renderVolumeButton.connect('clicked()', self.renderVolume)
    self.calcTubePositionButton.connect('clicked()', self.calcTubePosition)
    
    #self.maxButton.connect("clicked()", self.maxButtonClicked)
    #self.minButton.connect("clicked()", self.minButtonClicked)
   
  
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
  

  
  '''    
    def calcTubeData(self):
    #segmentationNode = self.segmentSelector.currentNode() 
    #slicer.modules.segmentations.logic().ExportVisibleSegmentsToLabelmapNode(segmentationNode, labelmapVolumeNode, self.volumeNode)
    segmentationNode = self.segmentSelector.currentNode() 
    
    #print(type(segmentationNode))
    
    segmentationNode.GetSegmentation().SetConversionParameter('Smoothing factor','0.0')
    segmentationNode.CreateClosedSurfaceRepresentation()
    polyData = segmentationNode.GetClosedSurfaceInternalRepresentation('Segment_2')
    
    points = polyData.GetPoints()
    pointsNum = points.GetNumberOfPoints()
    
    
    for point_id in range(pointsNum):
        pt=[]
        coords = points.GetPoint(point_id)
        f.write("{} {} {}\n".format(coords[0], coords[1], coords[2]))
    
    #rint(scals)  
    
    #print("Point number:", pointNum)
      

    
    
  def applySegment(self):
    # Input nodes
    #volumeNode = getNode('MRHead')
    segmentationNode = self.segmentSelector.currentNode()  
    self.volumeNode = self.volumeSelector.currentNode()

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
    maskedVolumeNode = slicer.modules.volumes.logic().CloneVolume(self.volumeNode, "Masked")
    slicer.util.updateVolumeFromArray(maskedVolumeNode, maskedVoxels)
    slicer.util.setSliceViewerLayers(maskedVolumeNode)
    
    self.volumeNode = self.volumeSelector.currentNode()
    
  def renderVolume(self):
    logic = slicer.modules.volumerendering.logic()
    self.volumeNode = self.volumeSelector.currentNode()
    self.volumePreset = logic.GetPresetByName(self.ctTypeSelector.currentText)
    
    if self.roiNode == None:
        self.roiNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLAnnotationROINode")
        self.roiNode.SetName("Masked render ROI")
        
    if self.displayNode!=None:
        slicer.mrmlScene.RemoveNode(self.displayNode)
        self.displayNode=None
    
    if self.volumeNode!=None:
        self.displayNode = logic.CreateVolumeRenderingDisplayNode()
        self.displayNode.SetAndObserveROINodeID(self.roiNode.GetID())
        self.displayNode.SetCroppingEnabled(1)
        self.displayNode.UnRegister(logic)
        print(self.displayNode)        
        slicer.mrmlScene.AddNode(self.displayNode)
        
        self.volumeNode.AddAndObserveDisplayNodeID(self.displayNode.GetID())
        logic.UpdateDisplayNodeFromVolumeNode(self.displayNode, self.volumeNode)
        self.propertyNode = self.displayNode.GetVolumePropertyNode()
        self.propertyNode.Copy(self.volumePreset)
        
        
        self.displayNode.SetVisibility(1)
        #self.propertyNode = self.displayNode.GetVolumePropertyNode()
    def applyROI(self):
      self.displayNode.SetAndObserveROINodeID(self.roiNode.GetID())
      self.displayNode.CroppingEnabled = 1 
  
  def createVolume(self):
    nodeName = "MyNewVolume"
    imageSize = [512, 512, 512]
    voxelType=vtk.VTK_UNSIGNED_CHAR
    imageOrigin = [0.0, 0.0, 0.0]
    imageSpacing = [1.0, 1.0, 1.0]
    imageDirections = [[1,0,0], [0,1,0], [0,0,1]]
    fillVoxelValue = 1.0

    # Create an empty image volume, filled with fillVoxelValue
    self.imageData = vtk.vtkImageData()
    self.imageData.SetDimensions(imageSize)
    self.imageData.AllocateScalars(voxelType, 1)
    self.thresholder = vtk.vtkImageThreshold()
    self.thresholder.SetInputData(self.imageData)
    self.thresholder.SetInValue(fillVoxelValue)
    self.thresholder.SetOutValue(fillVoxelValue)
    self.thresholder.Update()
    # Create volume node
    self.volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", nodeName)
    self.volumeNode.SetOrigin(imageOrigin)
    self.volumeNode.SetSpacing(imageSpacing)
    self.volumeNode.SetIJKToRASDirections(imageDirections)
    self.volumeNode.SetAndObserveImageData(self.thresholder.GetOutput())
    self.volumeNode.CreateDefaultDisplayNodes()
    self.volumeNode.CreateDefaultStorageNode()
    
    ijkToRas = vtk.vtkMatrix4x4()
    self.volumeNode.GetIJKToRASMatrix(ijkToRas)
    imageData=self.volumeNode.GetImageData()
    extent = imageData.GetExtent()
    for k in range(extent[4], extent[5]+1):
      print("k=",k)
      for j in range(extent[2], extent[3]+1):
        print("j=",j)
        for i in range(extent[0], extent[1]+1):
          position_Ijk=[i, j, k, 1]
          position_Ras=ijkToRas.MultiplyPoint(position_Ijk)
          r=position_Ras[0]
          a=position_Ras[1]
          s=position_Ras[2]      
          functionValue=(r-10)*(r-10)+(a+15)*(a+15)+s*s
          imageData.SetScalarComponentFromDouble(i,j,k,0,functionValue)
    imageData.SetScalarComponentFromFloat(distortionVectorPosition_Ijk[0], distortionVectorPosition_Ijk[1], distortionVectorPosition_Ijk[2], 0, fillValue)
    imageData.Modified()
     
    
  def setMinMaxVolumeSettings(self, min_val, max_val, thresh):
    of = vtk.vtkPiecewiseFunction()
    of.AddPoint(-10000, thresh)
    of.AddPoint(min_val, thresh)
    of.AddPoint(min_val, 0.0)
    of.AddPoint(max_val, 0.0)
    of.AddPoint(max_val, thresh)
    of.AddPoint(10000, thresh)
    self.propertyNode = self.displayNode.GetVolumePropertyNode()
    self.propertyNode.SetScalarOpacity(of)
        
  def minButtonClicked(self):
    min_v = int(self.lineEditMin.text)
    max_v = int(self.lineEditMax.text)
    print("min = ",min_v, " max=", max_v)
    self.setMinMaxVolumeSettings(min_v, max_v, 1.0)
    
    #"Segment_1"

    self.markupSelectorLabel = qt.QLabel()
    self.markupSelectorLabel.setText( "Markup list: " )
    self.markupSelector = slicer.qMRMLNodeComboBox()
    self.markupSelector.nodeTypes = ( "vtkMRMLMarkupsFiducialNode", "" )
    self.markupSelector.noneEnabled = False
    self.markupSelector.selectNodeUponCreation = True
    self.markupSelector.setMRMLScene( slicer.mrmlScene )
    self.markupSelector.setToolTip( "Pick the markup list to be filled" )
    layout.addRow(self.markupSelectorLabel, self.markupSelector) 
      
        
    # Apply button
    self.computeButton = qt.QPushButton("Compute")
    self.computeButton.toolTip = "Compute information for the selected markup"    
    layout.addWidget(self.computeButton)
    self.UpdatecomputeButtonState()
    
    # Results
    self.totalDistanceLabel = qt.QLabel()
    self.totalDistanceLabel.setText( "Total distance between fiducials (mm): " )
    self.totalDistanceValue = qt.QLabel()
    layout.addRow(self.totalDistanceLabel, self.totalDistanceValue)    

    # connections
    self.computeButton.connect('clicked()', self.onCompute)
    self.markupSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.onMarkupSelect)    
    ''' 
  
  '''
  def UpdatecomputeButtonState(self):
    if not self.markupSelector.currentNode() :
      self.computeButton.enabled = False
    else:
      self.computeButton.enabled = True      
    
  def onMarkupSelect(self, node):
    self.UpdatecomputeButtonState()

  def onCompute(self):
    slicer.app.processEvents()
    self.logic = MarkupsInfoLogic(self.markupSelector.currentNode())
    self.totalDistanceValue.setText('%.2f'%self.logic.info['totalDistance'])
  '''

'''

#
# Logic
#
    
class MarkupsInfoLogic:
  """Implement the logic to compute markup info
  Nodes are passed in as arguments.
  Results are stored as 'info' instance variable.
  """
  
  def __init__(self, markupNode):    
    
    self.info={}
    
    # Compute total distance between fiducials        
    totalDist=0
    startPtCoords = [0.0, 0.0, 0.0]
    endPtCoords = [0.0, 0.0, 0.0]
    for fidIndex in range(markupNode.GetNumberOfFiducials()-1): 
        markupNode.GetNthFiducialPosition(fidIndex,startPtCoords)
        markupNode.GetNthFiducialPosition(fidIndex+1,endPtCoords)
        dist=math.sqrt((startPtCoords[0]-endPtCoords[0])**2+(startPtCoords[1]-endPtCoords[1])**2+(startPtCoords[2]-endPtCoords[2])**2)
        totalDist=totalDist+dist
    self.info['totalDistance']=totalDist
'''
  
