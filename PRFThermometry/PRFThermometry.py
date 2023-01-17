import os
import unittest
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import SimpleITK as sitk
import sitkUtils
import numpy
import math
import copy
from skimage.restoration import unwrap_phase

#
# PRFThermometry
#

class PRFThermometry(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "PRFThermometry" # TODO make this more human readable by adding spaces
    self.parent.categories = ["IGT"]
    self.parent.dependencies = []
    self.parent.contributors = ["Junichi Tokuda (Brigham and Women's Hospital)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This module was developed based on a template created by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# PRFThermometryWidget
#

class PRFThermometryWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):

    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    # --------------------------------------------------
    # I/O Area
    # --------------------------------------------------
    ioCollapsibleButton = ctk.ctkCollapsibleButton()
    ioCollapsibleButton.text = "I/O"
    self.layout.addWidget(ioCollapsibleButton)

    ioFormLayout = qt.QVBoxLayout(ioCollapsibleButton)
    
    ioCommonFormLayout = qt.QFormLayout()
    ioFormLayout.addLayout(ioCommonFormLayout)

    # --------------------
    # Single Frame
    #
    singleFrameGroupBox = ctk.ctkCollapsibleGroupBox()
    singleFrameGroupBox.title = "Single Frame"
    singleFrameGroupBox.collapsed = False

    ioFormLayout.addWidget(singleFrameGroupBox)
    singleFrameFormLayout = qt.QFormLayout(singleFrameGroupBox)
    
    #
    # input volume selector
    #
    self.baselinePhaseSelector = slicer.qMRMLNodeComboBox()
    self.baselinePhaseSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.baselinePhaseSelector.selectNodeUponCreation = True
    self.baselinePhaseSelector.addEnabled = False
    self.baselinePhaseSelector.removeEnabled = False
    self.baselinePhaseSelector.noneEnabled = True
    self.baselinePhaseSelector.showHidden = False
    self.baselinePhaseSelector.showChildNodeTypes = False
    self.baselinePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.baselinePhaseSelector.setToolTip( "Pick the baseline phase map" )
    singleFrameFormLayout.addRow("Baseline Phase Volume: ", self.baselinePhaseSelector)

    #
    # input volume selector
    #
    self.referencePhaseSelector = slicer.qMRMLNodeComboBox()
    self.referencePhaseSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.referencePhaseSelector.selectNodeUponCreation = True
    self.referencePhaseSelector.addEnabled = False
    self.referencePhaseSelector.removeEnabled = False
    self.referencePhaseSelector.noneEnabled = True
    self.referencePhaseSelector.showHidden = False
    self.referencePhaseSelector.showChildNodeTypes = False
    self.referencePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.referencePhaseSelector.setToolTip( "Select a reference phase volume." )
    singleFrameFormLayout.addRow("Reference Phase Volume: ", self.referencePhaseSelector)

    #
    # input volume selector
    #
    self.maskSelector = slicer.qMRMLNodeComboBox()
    self.maskSelector.nodeTypes = ( ("vtkMRMLLabelMapVolumeNode"), "" )
    self.maskSelector.selectNodeUponCreation = True
    self.maskSelector.addEnabled = False
    self.maskSelector.removeEnabled = False
    self.maskSelector.noneEnabled = True
    self.maskSelector.showHidden = False
    self.maskSelector.showChildNodeTypes = False
    self.maskSelector.setMRMLScene( slicer.mrmlScene )
    self.maskSelector.setToolTip( "Select a mask volume." )
    singleFrameFormLayout.addRow("Mask Volume: ", self.maskSelector)
    
    #
    # tempMap volume selector
    #
    self.tempMapSelector = slicer.qMRMLNodeComboBox()
    self.tempMapSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.tempMapSelector.selectNodeUponCreation = True
    self.tempMapSelector.addEnabled = True
    self.tempMapSelector.removeEnabled = True
    self.tempMapSelector.noneEnabled = True
    self.tempMapSelector.renameEnabled = True
    self.tempMapSelector.showHidden = False
    self.tempMapSelector.showChildNodeTypes = False
    self.tempMapSelector.setMRMLScene( slicer.mrmlScene )
    self.tempMapSelector.setToolTip( "Select an output temperature map." )
    singleFrameFormLayout.addRow("Output Temperature Map: ", self.tempMapSelector)

    #
    # Apply Button
    #
    self.applyButtonSingle = qt.QPushButton("Generate Single-Frame Temperature Map")
    self.applyButtonSingle.toolTip = "Run the algorithm for the single frame input."
    self.applyButtonSingle.enabled = False
    singleFrameFormLayout.addRow(self.applyButtonSingle)


    # --------------------
    # Multi-Frame
    #
    multiFrameGroupBox = ctk.ctkCollapsibleGroupBox()
    multiFrameGroupBox.title = "Multi Frame"
    multiFrameGroupBox.collapsed = True

    ioFormLayout.addWidget(multiFrameGroupBox)
    multiFrameFormLayout = qt.QFormLayout(multiFrameGroupBox)
    
    #
    # input volume selector
    #
    self.multiFrameReferencePhaseSelector = slicer.qMRMLNodeComboBox()
    self.multiFrameReferencePhaseSelector.nodeTypes = ( ("vtkMRMLSequenceNode"), "" )
    self.multiFrameReferencePhaseSelector.selectNodeUponCreation = True
    self.multiFrameReferencePhaseSelector.addEnabled = False
    self.multiFrameReferencePhaseSelector.removeEnabled = False
    self.multiFrameReferencePhaseSelector.noneEnabled = False
    self.multiFrameReferencePhaseSelector.showHidden = False
    self.multiFrameReferencePhaseSelector.showChildNodeTypes = False
    self.multiFrameReferencePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.multiFrameReferencePhaseSelector.setToolTip( "Select a sequence node that contains the reference phase maps" )
    multiFrameFormLayout.addRow("Reference phase sequence: ", self.multiFrameReferencePhaseSelector)

    #
    # tempMap volume selector
    #
    self.multiFrameTempMapSelector = slicer.qMRMLNodeComboBox()
    self.multiFrameTempMapSelector.nodeTypes = ( ("vtkMRMLSequenceNode"), "" )
    self.multiFrameTempMapSelector.selectNodeUponCreation = True
    self.multiFrameTempMapSelector.addEnabled = True
    self.multiFrameTempMapSelector.removeEnabled = True
    self.multiFrameTempMapSelector.noneEnabled = True
    self.multiFrameTempMapSelector.renameEnabled = True
    self.multiFrameTempMapSelector.showHidden = False
    self.multiFrameTempMapSelector.showChildNodeTypes = False
    self.multiFrameTempMapSelector.setMRMLScene( slicer.mrmlScene )
    self.multiFrameTempMapSelector.setToolTip( "Select an output sequence to store temperature maps." )
    multiFrameFormLayout.addRow("Output Temperature Map: ", self.multiFrameTempMapSelector)
 
    #
    # Apply Button
    #
    self.applyButtonMulti = qt.QPushButton("Generate Multi-Frame Temperature Map")
    self.applyButtonMulti.toolTip = "Run the algorithm for the multi-frame input."
    self.applyButtonMulti.enabled = False
    multiFrameFormLayout.addRow(self.applyButtonMulti)

    
    # --------------------
    # Simple Masking
    #
    maskingGroupBox = ctk.ctkCollapsibleGroupBox()
    maskingGroupBox.title = "Simple Masking"
    maskingGroupBox.collapsed = True

    ioFormLayout.addWidget(maskingGroupBox)
    maskingFormLayout = qt.QFormLayout(maskingGroupBox)
    
    self.simpleMaskingFlagCheckBox = qt.QCheckBox()
    self.simpleMaskingFlagCheckBox.checked = 0
    self.simpleMaskingFlagCheckBox.setToolTip("If checked, simple masking is ON.")
    maskingFormLayout.addRow("Use Simple Masking", self.simpleMaskingFlagCheckBox)
    
    self.radiusSpinBox = qt.QDoubleSpinBox()
    self.radiusSpinBox.objectName = 'radiusSpinBox'
    self.radiusSpinBox.setMaximum(10.0)
    self.radiusSpinBox.setMinimum(0.0)
    self.radiusSpinBox.setDecimals(4)
    self.radiusSpinBox.setValue(0.8)
    self.radiusSpinBox.setToolTip("Radius of the disk mask")
    maskingFormLayout.addRow("Radius: ", self.radiusSpinBox)


    # --------------------------------------------------
    # PRF Parameters Area
    # --------------------------------------------------
    #
    prfParametersCollapsibleButton = ctk.ctkCollapsibleButton()
    prfParametersCollapsibleButton.text = "PRF Parameters"
    self.layout.addWidget(prfParametersCollapsibleButton)

    prfParametersFormLayout = qt.QFormLayout(prfParametersCollapsibleButton)
    
    #
    # Temperature coefficient (alpha)
    #
    self.alphaSpinBox = qt.QDoubleSpinBox()
    self.alphaSpinBox.objectName = 'alphaSpinBox'
    self.alphaSpinBox.setMaximum(100.0)
    self.alphaSpinBox.setMinimum(-100.0)
    self.alphaSpinBox.setDecimals(8)
    self.alphaSpinBox.setValue(-0.01)
    self.alphaSpinBox.setToolTip("Temperature coefficient (ppm/deg C)")
    prfParametersFormLayout.addRow("alpha (ppm/deg C): ", self.alphaSpinBox)

    #
    # Gyromagnetic ratio (gamma)
    #
    self.gammaSpinBox = qt.QDoubleSpinBox()
    self.gammaSpinBox.objectName = 'gammaSpinBox'
    self.gammaSpinBox.setMaximum(100.0)
    self.gammaSpinBox.setMinimum(-100.0)
    self.gammaSpinBox.setDecimals(8)
    self.gammaSpinBox.setValue(42.576)
    self.gammaSpinBox.setToolTip("Gyromagnetic ratio / PI (rad / s T) (Default is proton )")
    prfParametersFormLayout.addRow("gamma (rad / s T): ", self.gammaSpinBox)

    #
    # Field strength (B0)
    #
    self.B0SpinBox = qt.QDoubleSpinBox()
    self.B0SpinBox.objectName = 'B0SpinBox'
    self.B0SpinBox.setMaximum(20.0)
    self.B0SpinBox.setMinimum(0.0)
    self.B0SpinBox.setDecimals(8)
    self.B0SpinBox.setValue(3.0)
    self.B0SpinBox.setToolTip("Static field strength (Tesla)")
    prfParametersFormLayout.addRow("B0 (Tesla): ", self.B0SpinBox)

    #
    # Echo time (TE)
    #
    self.TESpinBox = qt.QDoubleSpinBox()
    self.TESpinBox.objectName = 'TESpinBox'
    self.TESpinBox.setMaximum(10.0)
    self.TESpinBox.setMinimum(0.0)
    self.TESpinBox.setDecimals(8)
    self.TESpinBox.setValue(0.01)
    self.TESpinBox.setToolTip("Echo time (s)")
    prfParametersFormLayout.addRow("TE (s): ", self.TESpinBox)

    #
    # Body temperature (BT)
    #
    self.BTSpinBox = qt.QDoubleSpinBox()
    self.BTSpinBox.objectName = 'BTSpinBox'
    self.BTSpinBox.setMaximum(100.0)
    self.BTSpinBox.setMinimum(0.0)
    self.BTSpinBox.setDecimals(8)
    self.BTSpinBox.setValue(37.0)
    self.BTSpinBox.setToolTip("Body temperature (Deg C)")
    prfParametersFormLayout.addRow("BT (Deg C): ", self.BTSpinBox)

    
    # --------------------------------------------------
    # Parameters Area
    # --------------------------------------------------
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # Check box to apply phase unwrapping
    #
    self.phaseUnwrappingFlagCheckBox = qt.QCheckBox()
    self.phaseUnwrappingFlagCheckBox.checked = 0
    self.phaseUnwrappingFlagCheckBox.setToolTip("If checked, use phase unwrapping on the raw input images.")
    parametersFormLayout.addRow("Phase unwrapping on raw input", self.phaseUnwrappingFlagCheckBox)
    
    self.phaseUnwrappingPostFlagCheckBox = qt.QCheckBox()
    self.phaseUnwrappingPostFlagCheckBox.checked = 1
    self.phaseUnwrappingPostFlagCheckBox.setToolTip("If checked, use phase unwrapping after computing the phase shift.")
    parametersFormLayout.addRow("Phase unwrapping after subtraction", self.phaseUnwrappingPostFlagCheckBox)

    self.complexFlagCheckBox = qt.QCheckBox()
    self.complexFlagCheckBox.checked = 1
    self.complexFlagCheckBox.setToolTip("If checked, use complex values to subtract phase.")
    parametersFormLayout.addRow("Use complex values", self.complexFlagCheckBox)

    #
    # Phase range (when "Use complex values" is ON)
    #
    self.phaseRangeSpinBox = qt.QDoubleSpinBox()
    self.phaseRangeSpinBox.objectName = 'alphaSpinBox'
    self.phaseRangeSpinBox.setMaximum(180.0)
    self.phaseRangeSpinBox.setMinimum(-180.0)
    self.phaseRangeSpinBox.setDecimals(2)
    self.phaseRangeSpinBox.setValue(30.0)
    self.phaseRangeSpinBox.setToolTip("Specify the phase range in degrees. The normal phase range for temperature mapping is [-360 deg, 0 deg]. If '30' is speicified, the range is shifted by 30 degrees (i.e., [-330deg, 30deg]). The shift will allow visualizing negative temperature changes even when complex values are used to calcluate the phase shift.")
    parametersFormLayout.addRow("Phase range shift (deg): ", self.phaseRangeSpinBox)
    
    #
    # Check box to use threshold
    #
    self.useThresholdFlagCheckBox = qt.QCheckBox()
    self.useThresholdFlagCheckBox.checked = 1
    self.useThresholdFlagCheckBox.setToolTip("If checked, apply the threshold to limit the pixel value ranges.")
    parametersFormLayout.addRow("Use Threshold", self.useThresholdFlagCheckBox)
    
    #
    # Temperature color range 
    #
    self.scaleRangeMaxSpinBox = qt.QDoubleSpinBox()
    self.scaleRangeMaxSpinBox.objectName = 'scaleRangeMaxSpinBox'
    self.scaleRangeMaxSpinBox.setMaximum(200.0)
    self.scaleRangeMaxSpinBox.setMinimum(0.0)
    self.scaleRangeMaxSpinBox.setDecimals(4)
    self.scaleRangeMaxSpinBox.setValue(85.0)
    self.scaleRangeMaxSpinBox.setToolTip("Maximum value for the temperature color scale.")
    parametersFormLayout.addRow("Max. scale (deg C): ", self.scaleRangeMaxSpinBox)

    #
    # Lower threshold - We set threshold value to limit the range of intensity 
    #
    self.scaleRangeMinSpinBox = qt.QDoubleSpinBox()
    self.scaleRangeMinSpinBox.objectName = 'scaleRangeMinSpinBox'
    self.scaleRangeMinSpinBox.setMaximum(200.0)
    self.scaleRangeMinSpinBox.setMinimum(0.0)
    self.scaleRangeMinSpinBox.setDecimals(4)
    self.scaleRangeMinSpinBox.setValue(35.0)
    self.scaleRangeMinSpinBox.setToolTip("Minimum value for the temperature color scale")
    parametersFormLayout.addRow("Min. scale (deg C): ", self.scaleRangeMinSpinBox)

    #
    # Upper threshold - We set threshold value to limit the range of intensity 
    #
    self.upperThresholdSpinBox = qt.QDoubleSpinBox()
    self.upperThresholdSpinBox.objectName = 'upperThresholdSpinBox'
    self.upperThresholdSpinBox.setMaximum(1000000.0)
    self.upperThresholdSpinBox.setMinimum(-1000000.0)
    self.upperThresholdSpinBox.setDecimals(6)
    self.upperThresholdSpinBox.setValue(1000.0)
    self.upperThresholdSpinBox.setToolTip("Upper threshold for the output")
    parametersFormLayout.addRow("Upper Threshold (deg): ", self.upperThresholdSpinBox)

    #
    # Lower threshold - We set threshold value to limit the range of intensity 
    #
    self.lowerThresholdSpinBox = qt.QDoubleSpinBox()
    self.lowerThresholdSpinBox.objectName = 'lowerThresholdSpinBox'
    self.lowerThresholdSpinBox.setMaximum(1000000.0)
    self.lowerThresholdSpinBox.setMinimum(-1000000.0)
    self.lowerThresholdSpinBox.setDecimals(6)
    self.lowerThresholdSpinBox.setValue(-1000.0)
    self.lowerThresholdSpinBox.setToolTip("Lower threshold for the output")
    parametersFormLayout.addRow("Lower Threshold (deg): ", self.lowerThresholdSpinBox)



    
    #
    # Check for automatic update
    #
    self.autoUpdateCheckBox = qt.QCheckBox()
    self.autoUpdateCheckBox.checked = 1
    self.autoUpdateCheckBox.setToolTip("Automatic Update")
    parametersFormLayout.addRow("Automatic Update", self.autoUpdateCheckBox)

    # connections
    
    self.applyButtonSingle.connect('clicked(bool)', self.onApplyButtonSingle)
    self.applyButtonMulti.connect("clicked(bool)", self.onApplyButtonMulti)
        
    self.baselinePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelectSingle)
    self.referencePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelectSingle)
    self.tempMapSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelectSingle)
    self.multiFrameReferencePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelectMulti)
    self.multiFrameTempMapSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelectMulti)

    self.complexFlagCheckBox.connect('toggled(bool)', self.onComplexFlag)
    
    self.useThresholdFlagCheckBox.connect('toggled(bool)', self.onUseThreshold)
    self.autoUpdateCheckBox.connect('toggled(bool)', self.onAutoUpdate)

    self.simpleMaskingFlagCheckBox.connect('toggled(bool)', self.onUseSimpleMask)
    
    
    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelectSingle()

    self.tag = None

    
  def cleanup(self):
    pass

  
  def onSelectSingle(self):
    self.applyButtonSingle.enabled = self.baselinePhaseSelector.currentNode() and self.baselinePhaseSelector.currentNode() and self.tempMapSelector.currentNode()

    
  def onSelectMulti(self):
    self.applyButtonMulti.enabled = self.multiFrameReferencePhaseSelector.currentNode() and self.multiFrameTempMapSelector.currentNode()


  def onComplexFlag(self):
    if self.complexFlagCheckBox.checked == True:
      self.phaseRangeSpinBox.enabled = True
    else:
      self.phaseRangeSpinBox.enabled = False
    
  def onUseThreshold(self):
    if self.useThresholdFlagCheckBox.checked == True:
      self.lowerThresholdSpinBox.enabled = True;      
      self.upperThresholdSpinBox.enabled = True;      
    else:
      self.lowerThresholdSpinBox.enabled = False;      
      self.upperThresholdSpinBox.enabled = False;      

  def onUseSimpleMask(self):
    
    if self.simpleMaskingFlagCheckBox.checked == 1:
      self.maskSelector.enabled = False
    else:
      self.maskSelector.enabled = True
      
      
  def onAutoUpdate(self):
    if self.autoUpdateCheckBox.checked == True:
      if self.referencePhaseSelector.currentNode():
        self.referencePhaseSelector.enabled = False
        refNode = self.referencePhaseSelector.currentNode()
        self.tag = refNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.onModelRefImageModifiedEvent)
      else: # Cannot set autoupdate 
        self.autoUpdateCheckBox.checked = Falase
    else:
      if self.tag:
        if self.referencePhaseSelector.currentNode():
          refNode = self.referencePhaseSelector.currentNode()
          refNode.RemoveObserver(self.tag)
        self.referencePhaseSelector.enabled = True

        
  def onModelRefImageModifiedEvent(self, caller, event):
    self.onApplyButtonSingle()

    
  def onApplyButtonSingle(self):
    logic = PRFThermometryLogic()

    param = {}
    param['usePhaseUnwrapping']       = self.phaseUnwrappingFlagCheckBox.checked
    param['usePhaseUnwrappingPost']   = self.phaseUnwrappingPostFlagCheckBox.checked
    param['useComplex']               = self.complexFlagCheckBox.checked
    param['phaseRangeShiftDeg']       = self.phaseRangeSpinBox.value
    param['baselinePhaseVolumeNode']  = self.baselinePhaseSelector.currentNode()
    param['referencePhaseVolumeNode'] = self.referencePhaseSelector.currentNode()
    param['maskVolumeNode']           = self.maskSelector.currentNode()
    param['tempMapVolumeNode']        = self.tempMapSelector.currentNode()
    param['alpha']                    = self.alphaSpinBox.value
    param['gamma']                    = self.gammaSpinBox.value
    param['B0']                       = self.B0SpinBox.value
    param['TE']                       = self.TESpinBox.value
    param['BT']                       = self.BTSpinBox.value
    param['colorScaleMax']            = self.scaleRangeMaxSpinBox.value
    param['colorScaleMin']            = self.scaleRangeMinSpinBox.value
    
    if self.useThresholdFlagCheckBox.checked == True:
      param['upperThreshold']         = self.upperThresholdSpinBox.value
      param['lowerThreshold']         = self.lowerThresholdSpinBox.value
    else:
      param['upperThreshold']         = False
      param['lowerThreshold']         = False

    if self.simpleMaskingFlagCheckBox.checked == 1:
      param['simpleMask']        = 'disk'
      param['simpleMask.radius'] = self.radiusSpinBox.value
    else:
      param['simpleMask']        = None

      
    logic.runSingleFrame(param)


  def onApplyButtonMulti(self):
    logic = PRFThermometryLogic()

    param = {}
    param['baselinePhaseVolumeNode']  = self.baselinePhaseSelector.currentNode()
    param['maskVolumeNode']           = self.maskSelector.currentNode()
    param['referencePhaseSequenceNode'] = self.multiFrameReferencePhaseSelector.currentNode()
    param['tempMapSequenceNode']        = self.multiFrameTempMapSelector.currentNode()
    param['usePhaseUnwrapping']       = self.phaseUnwrappingFlagCheckBox.checked
    param['usePhaseUnwrappingPost']   = self.phaseUnwrappingPostFlagCheckBox.checked
    param['useComplex']               = self.complexFlagCheckBox.checked
    param['phaseRangeShiftDeg']       = self.phaseRangeSpinBox.value
    #param['truePhasePointNode']       = self.truePhasePointSelector.currentNode()
    param['alpha']                    = self.alphaSpinBox.value
    param['gamma']                    = self.gammaSpinBox.value
    param['B0']                       = self.B0SpinBox.value
    param['TE']                       = self.TESpinBox.value
    param['BT']                       = self.BTSpinBox.value
    param['colorScaleMax']            = self.scaleRangeMaxSpinBox.value
    param['colorScaleMin']            = self.scaleRangeMinSpinBox.value
    
    if self.useThresholdFlagCheckBox.checked == True:
      param['upperThreshold']         = self.upperThresholdSpinBox.value
      param['lowerThreshold']         = self.lowerThresholdSpinBox.value
    else:
      param['upperThreshold']         = False
      param['lowerThreshold']         = False

    logic.runMultiFrame(param)


    
  def onReload(self, moduleName="PRFThermometry"):
    # Generic reload method for any scripted module.
    # ModuleWizard will subsitute correct default moduleName.

    globals()[moduleName] = slicer.util.reloadScriptedModule(moduleName)


#
# PRFThermometryLogic
#

class PRFThermometryLogic(ScriptedLoadableModuleLogic):

  def isValidInputOutputData(self, baselinePhaseVolumeNode, referencePhaseVolumeNode):
    """Validates if the output is not the same as input
    """
    if not baselinePhaseVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node for ParamA image defined')
      return False
    if not referencePhaseVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node for ParamB image defined')
      return False
    return True


  def generateDiskMask(self, refImage, center=[0.5,0.5,0.5], radius=0.5):

    refImageArray = sitk.GetArrayFromImage(refImage)
    dims = refImageArray.shape
    maxDim = numpy.double(numpy.max(dims))
    r = maxDim * radius
    
    cx = dims[0] * center[0]
    cy = dims[1] * center[1]
    cz = dims[2] * center[2]
    y, x, z = numpy.ogrid[-cx:dims[0]-cx, -cy:dims[1]-cy, -cz:dims[2]-cz]
    mask = x*x + y*y + z*z <= r*r
    voxels = numpy.zeros((dims[0], dims[1], dims[2]))
    voxels[mask] = 1
    
    maskImage = sitk.GetImageFromArray(voxels)
    maskImage.CopyInformation(refImage) 
    
    return maskImage

  
  def runSingleFrame(self, param): 
    """
    Run the actual algorithm
    """
    usePhaseUnwrapping       = param['usePhaseUnwrapping']
    usePhaseUnwrappingPost   = param['usePhaseUnwrappingPost']
    useComplex               = param['useComplex']
    phaseRangeShiftDeg       = param['phaseRangeShiftDeg']
    baselinePhaseVolumeNode  = param['baselinePhaseVolumeNode']
    referencePhaseVolumeNode = param['referencePhaseVolumeNode']
    maskVolumeNode           = param['maskVolumeNode']
    tempMapVolumeNode        = param['tempMapVolumeNode']
    #truePhasePointNode       = param['truePhasePointNode']
    alpha                    = param['alpha']
    gamma                    = param['gamma']
    B0                       = param['B0']
    TE                       = param['TE']
    BT                       = param['BT']
    colorScaleMax            = param['colorScaleMax']
    colorScaleMin            = param['colorScaleMin']
    upperThreshold           = param['upperThreshold']
    lowerThreshold           = param['lowerThreshold']

    # Convert the phase range shift from degree to radian
    phaseRangeShift = numpy.pi * phaseRangeShiftDeg/180.0
    
    if not self.isValidInputOutputData(baselinePhaseVolumeNode, referencePhaseVolumeNode):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    print(baselinePhaseVolumeNode.GetName())
    print(referencePhaseVolumeNode.GetName())

    imageBaseline = None
    imageReference = None
    
    # Check the scalar type (Siemens SRC sends image data in 'short' instead of 'unsigned short')
    baselineImageData = baselinePhaseVolumeNode.GetImageData()
    scalarType = ''
    if baselineImageData != None:
      scalarType = baselineImageData.GetScalarTypeAsString()
      
    imageBaseline  = sitk.Cast(sitkUtils.PullVolumeFromSlicer(baselinePhaseVolumeNode), sitk.sitkFloat64)
    imageReference  = sitk.Cast(sitkUtils.PullVolumeFromSlicer(referencePhaseVolumeNode), sitk.sitkFloat64)
    
    mask = None
    
    if param['simpleMask'] == 'disk':
      simpleMaskRadius         = param['simpleMask.radius'] 
      mask = self.generateDiskMask(imageBaseline, radius=simpleMaskRadius)
      imageBaseline = imageBaseline * mask
      imageReference = imageReference * mask
    elif param['maskVolumeNode']:
      mask = sitk.Cast(sitkUtils.PullVolumeFromSlicer(maskVolumeNode), sitk.sitkFloat64)
      imageBaseline = imageBaseline * mask
      imageReference = imageReference * mask
    
    if tempMapVolumeNode:

      self.phaseDiff = None
      self.phaseDrift = None
      
      if scalarType == 'unsigned short':
        print('imageBaseline*numpy.pi/2048.0 - numpy.pi')
        imageBaselinePhase = imageBaseline*numpy.pi/2048.0 - numpy.pi
        imageReferencePhase = imageReference*numpy.pi/2048.0 - numpy.pi
      else:
        print('imageBaseline*numpy.pi/4096.0')
        imageBaselinePhase = imageBaseline*numpy.pi/4096.0
        imageReferencePhase = imageReference*numpy.pi/4096.0

      # Phase unwrapping on the raw input images
      if usePhaseUnwrapping == True:
        print('usePhaseUnwrapping')
        imageBaselinePhase  = self.unwrap(imageBaselinePhase)
        imageReferencePhase = self.unwrap(imageReferencePhase)

        # Match the phases
        # Phase unwrapping often ends up shifting the entire phase map by pi*N. Try shifting the reference phase map
        # by pi*N where N = -2, -1, ... ,2 and check the difference between the baseline and reference images.
        stat = sitk.StatisticsImageFilter()
        meanDiffList = numpy.array([])
        nList = list(range(-4,4))
        for n in nList:
          phaseShift = numpy.pi * n 
          diffImage = imageReferencePhase + phaseShift - imageBaselinePhase
          stat.Execute(diffImage)
          meanDiffList = numpy.append(meanDiffList, numpy.abs(stat.GetMean()))
        
        minIndex = numpy.argmin(meanDiffList)
        print('minIndex = ' + str(minIndex))
        phaseShift = numpy.pi * nList[minIndex]
        imageReferencePhase = imageReferencePhase + phaseShift

      # Phase shift can be determined by either subtracting the phase values or
      # calculating the rotation in the complex space
      if useComplex == True:
        # Convert the phase images to complex images (as numpy arrays)
        arrayBaseline = sitk.GetArrayFromImage(imageBaselinePhase)
        arrayBaselineComplex = numpy.cos(arrayBaseline) + numpy.sin(arrayBaseline) * 1.0j
        arrayReference = sitk.GetArrayFromImage(imageReferencePhase)
        arrayReferenceComplex = numpy.cos(arrayReference) + numpy.sin(arrayReference) * 1.0j
        
        arrayPhaseDiffComplex = arrayReferenceComplex / arrayBaselineComplex
        arrayPhaseDiff = numpy.angle(arrayPhaseDiffComplex)
        
        # Change the range from [-pi, pi] to [-2pi+numpy.pi/4, numpy.pi/4] (allow some temperature decrease)
        arrayPhaseDiff[arrayPhaseDiff>phaseRangeShift] -= 2*numpy.pi
        
        self.phaseDiff = sitk.GetImageFromArray(arrayPhaseDiff)
        self.phaseDiff.SetOrigin(imageBaseline.GetOrigin())
        self.phaseDiff.SetSpacing(imageBaseline.GetSpacing())
        self.phaseDiff.SetDirection(imageBaseline.GetDirection())
        
      else:
        self.phaseDiff = imageReferencePhase - imageBaselinePhase        
      
      if usePhaseUnwrappingPost == True:
        print('usePhaseUnwrappingPost')
        self.phaseDiff  = self.unwrap(self.phaseDiff)
        
      self.phaseDiff = self.phaseDiff
      print("(alpha, gamma, B0, TE, TE) = (%f, %f, %f, %f, %f)" % (alpha, gamma, B0, TE, BT))
      imageTemp = self.phaseDiff / (alpha * 2.0 * numpy.pi * gamma * B0 * TE) + BT
        
      if upperThreshold or lowerThreshold:
        imageTempThreshold = sitk.Threshold(imageTemp, lowerThreshold, upperThreshold, 0.0)
        sitkUtils.PushVolumeToSlicer(imageTempThreshold, tempMapVolumeNode.GetName(), 0, True)
      else:
        sitkUtils.PushVolumeToSlicer(imageTemp, tempMapVolumeNode.GetName(), 0, True)

      dnode = tempMapVolumeNode.GetDisplayNode()
      if dnode == None:
        dnode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLScalarVolumeDisplayNode')
        slicer.mrmlScene.AddNode(dnode)
        tempMapVolumeNode.SetAndObserveDisplayNodeID(dnode.GetID())
        
      dnode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileColdToHotRainbow.txt')
      dnode.SetWindowLevelLocked(0)
      dnode.SetAutoWindowLevel(0)
      dnode.SetWindowLevelMinMax(colorScaleMin, colorScaleMax)

      colorLegendDisplayNode = slicer.modules.colors.logic().AddDefaultColorLegendDisplayNode(tempMapVolumeNode)
      colorLegendDisplayNode.VisibilityOn()

    logging.info('Processing completed')

    return True

  
  def runMultiFrame(self, param):
    
    # Run the algorithm for each volume frame under the given sequence node.
    # If a baseline volume is not specified, the first frame is used as the baseline for the temperature calculation.
    # Note that this function temporarily copies the baseline and reference nodes from the sequence node
    # (which has its own scene) to the main Slicer scene before calling runSingleFrame(), and remove them
    # once the temperature map is calculated.
    
    refSeqNode     = param['referencePhaseSequenceNode']
    tempMapSeqNode = param['tempMapSequenceNode']

    nVolumes = refSeqNode.GetNumberOfDataNodes()
    singleParam = copy.copy(param)

    unit = refSeqNode.GetIndexUnit()
    prefix = '%s_TempMap_' % refSeqNode.GetName()

    # Set up the output sequence node
    tempMapSeqNode.SetIndexType(refSeqNode.GetIndexType())
    tempMapSeqNode.SetIndexName(refSeqNode.GetIndexName())
    tempMapSeqNode.SetIndexUnit(unit)
    
    # If 'baselinePhaseVolumeNode' is None, use the first image as a baseline
    if param['baselinePhaseVolumeNode'] == None:
      baselinePhaseVolumeNode = refSeqNode.GetNthDataNode(0)
      # Copy the node to the Slicer scene
      copiedBaselinePhaseVolumeNode = slicer.mrmlScene.CopyNode(baselinePhaseVolumeNode)
      singleParam['baselinePhaseVolumeNode']  = copiedBaselinePhaseVolumeNode
    
    # Create a tempMapVolumeNode.
    # We only create a single output TempMap node. The TempMap node will be deep-copied when
    # the node is added to the TempMap sequence node.
    tempMapNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLScalarVolumeNode')
    tempMapNode.SetName(prefix+'Temp')
    slicer.mrmlScene.AddNode(tempMapNode)
    singleParam['tempMapVolumeNode'] = tempMapNode
    
    for i in range(nVolumes):
      print('Processing image # %d / %d' % ((i+1), nVolumes))
      phaseVolumeNode = refSeqNode.GetNthDataNode(i)
      copiedPhaseVolumeNode = slicer.mrmlScene.CopyNode(phaseVolumeNode)
      indexValue = refSeqNode.GetNthIndexValue(i)
      singleParam['referencePhaseVolumeNode'] = copiedPhaseVolumeNode
      self.runSingleFrame(singleParam)

      tempMapSeqNode.SetDataNodeAtValue(tempMapNode, indexValue)
      n = tempMapSeqNode.GetItemNumberFromIndexValue(indexValue)
      # Give a new new to the copied TempMap under the sequence.
      dnode = tempMapSeqNode.GetNthDataNode(n)
      dnode.SetName('%s%s%s' % (prefix, indexValue, unit))

      slicer.mrmlScene.RemoveNode(copiedPhaseVolumeNode)

    if param['baselinePhaseVolumeNode'] == None:
      slicer.mrmlScene.RemoveNode(copiedBaselinePhaseVolumeNode)
    slicer.mrmlScene.RemoveNode(tempMapNode)
    
    colorScaleMax            = param['colorScaleMax']
    colorScaleMin            = param['colorScaleMin']
    
    self.setProxyNode(tempMapSeqNode, colorScaleMin, colorScaleMax)


  def setProxyNode(self, sequenceNode, scaleMin, scaleMax):

    # Find the first sequence browser node
    sbNode = None
    col = slicer.mrmlScene.GetNodesByClass('vtkMRMLSequenceBrowserNode')
    if col.GetNumberOfItems() > 0:
      # Always use the first node
      sbNode = col.GetItemAsObject(0)
    else:
      sbNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLSequenceBrowserNode')
      slicer.mrmlScene.AddNode(sbNode)

    # Add the sequence node as a synchronized sequence node
    postfix = sbNode.AddSynchronizedSequenceNode(sequenceNode)

    pNode = sbNode.GetProxyNode(sequenceNode)
    if pNode == None:
      # Create a new node
      pNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLScalarVolumeNode')
      pNode.SetName(sequenceNode.GetName())
      slicer.mrmlScene.AddNode(pNode)

    dNode = pNode.GetDisplayNode()
    if dNode == None:
      dNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLScalarVolumeDisplayNode')
      slicer.mrmlScene.AddNode(dNode)
      pNode.SetAndObserveDisplayNodeID(dNode.GetID())
      
    dNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileColdToHotRainbow.txt')
    dNode.SetWindowLevelLocked(0)
    dNode.SetAutoWindowLevel(0)
    dNode.SetWindowLevelMinMax(scaleMin, scaleMax)
    
  
  def unwrap(self, imagePhase):
    imagePhaseNP = sitk.GetArrayFromImage(imagePhase)
    imageUnwrappedNP = unwrap_phase(imagePhaseNP)
    imageUnwrapped = sitk.GetImageFromArray(imageUnwrappedNP)
    imageUnwrapped.SetOrigin(imagePhase.GetOrigin())
    imageUnwrapped.SetSpacing(imagePhase.GetSpacing())
    imageUnwrapped.SetDirection(imagePhase.GetDirection())
    
    return imageUnwrapped
    

class PRFThermometryTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_PRFThermometry1()

  def test_PRFThermometry1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    pass

    #self.delayDisplay("Starting the test")
    ##
    ## first, get some data
    ##
    #import urllib
    #downloads = (
    #    ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
    #    )
    #
    #for url,name,loader in downloads:
    #  filePath = slicer.app.temporaryPath + '/' + name
    #  if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
    #    logging.info('Requesting download %s from %s...\n' % (name, url))
    #    urllib.urlretrieve(url, filePath)
    #  if loader:
    #    logging.info('Loading %s...' % (name,))
    #    loader(filePath)
    #self.delayDisplay('Finished with download and loading')
    #
    #volumeNode = slicer.util.getNode(pattern="FA")
    #logic = PRFThermometryLogic()
    #self.assertTrue( logic.hasImageData(volumeNode) )
    #self.delayDisplay('Test passed!')
