import math
from __main__ import qt, slicer

#
# MarkupsInfo module
#

class MarkupsInfo:
  def __init__(self, parent):
    import string
    parent.title = "Markups info"
    parent.categories = ["Informatics"]
    parent.contributors = ["Andras Lasso (PerkLab)"]
    parent.helpText = string.Template("""
Use this module to get information about markups. The only metric it computes now is the total distance between the markup points in the selected list.
    """).substitute({ 'a':parent.slicerWikiUrl, 'b':slicer.app.majorVersion, 'c':slicer.app.minorVersion })
    parent.acknowledgementText = """
    Supported by SparKit and the Slicer Community. See http://www.slicerrt.org for details.
    """
    self.parent = parent

#
# Widget
#

class MarkupsInfoWidget:

  def __init__(self, parent=None):
    self.parent = parent
    self.logic = None
    
  def setup(self):

    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )
    
    # Markup selector
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
  
