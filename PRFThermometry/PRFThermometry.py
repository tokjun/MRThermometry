import os
import unittest
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import SimpleITK as sitk
import sitkUtils
import numpy

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

    #--------------------------------------------------
    # For debugging
    #
    # Reload and Test area
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload && Test"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    reloadCollapsibleButton.collapsed = True
    
    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "PRFThermometry Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)
    #
    #--------------------------------------------------

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # Check box to correct noise
    #
    self.rowPhaseImageFlagCheckBox = qt.QCheckBox()
    self.rowPhaseImageFlagCheckBox.checked = 1
    self.rowPhaseImageFlagCheckBox.setToolTip("If checked, use original phase images as inputs")
    parametersFormLayout.addRow("Use raw phase images", self.rowPhaseImageFlagCheckBox)

    #
    # input volume selector
    #
    self.baselinePhaseSelector = slicer.qMRMLNodeComboBox()
    self.baselinePhaseSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.baselinePhaseSelector.selectNodeUponCreation = True
    self.baselinePhaseSelector.addEnabled = False
    self.baselinePhaseSelector.removeEnabled = False
    self.baselinePhaseSelector.noneEnabled = False
    self.baselinePhaseSelector.showHidden = False
    self.baselinePhaseSelector.showChildNodeTypes = False
    self.baselinePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.baselinePhaseSelector.setToolTip( "Pick the baseline phase map" )
    parametersFormLayout.addRow("Baseline phase image: ", self.baselinePhaseSelector)

    #
    # input volume selector
    #
    self.referencePhaseSelector = slicer.qMRMLNodeComboBox()
    self.referencePhaseSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.referencePhaseSelector.selectNodeUponCreation = True
    self.referencePhaseSelector.addEnabled = False
    self.referencePhaseSelector.removeEnabled = False
    self.referencePhaseSelector.noneEnabled = False
    self.referencePhaseSelector.showHidden = False
    self.referencePhaseSelector.showChildNodeTypes = False
    self.referencePhaseSelector.setMRMLScene( slicer.mrmlScene )
    self.referencePhaseSelector.setToolTip( "Pick the reference phase map" )
    parametersFormLayout.addRow("Reference phase image: ", self.referencePhaseSelector)


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
    self.tempMapSelector.setToolTip( "Pick the output temperature map." )
    parametersFormLayout.addRow("Output Temperature Map: ", self.tempMapSelector)

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
    parametersFormLayout.addRow("alpha: ", self.alphaSpinBox)

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
    parametersFormLayout.addRow("gamma: ", self.gammaSpinBox)

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
    parametersFormLayout.addRow("B0: ", self.B0SpinBox)

    #
    # Echo time (TE)
    #
    self.TESpinBox = qt.QDoubleSpinBox()
    self.TESpinBox.objectName = 'TESpinBox'
    self.TESpinBox.setMaximum(10.0)
    self.TESpinBox.setMinimum(0.0)
    self.TESpinBox.setDecimals(8)
    self.TESpinBox.setValue(0.01)
    self.TESpinBox.setToolTip("Echo time (ms)")
    parametersFormLayout.addRow("TE: ", self.TESpinBox)

    #
    # Body temperature (BT)
    #
    self.BTSpinBox = qt.QDoubleSpinBox()
    self.BTSpinBox.objectName = 'BTSpinBox'
    self.BTSpinBox.setMaximum(100.0)
    self.BTSpinBox.setMinimum(0.0)
    self.BTSpinBox.setDecimals(8)
    self.BTSpinBox.setValue(37.0)
    self.BTSpinBox.setToolTip("Echo time (Deg C)")
    parametersFormLayout.addRow("BT: ", self.BTSpinBox)

    #
    # Check box to use threshold
    #
    self.useThresholdFlagCheckBox = qt.QCheckBox()
    self.useThresholdFlagCheckBox.checked = 1
    self.useThresholdFlagCheckBox.setToolTip("If checked, apply the threshold to limit the pixel value ranges.")
    parametersFormLayout.addRow("Use Threshold", self.useThresholdFlagCheckBox)
    
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
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.rowPhaseImageFlagCheckBox.connect('toggled(bool)', self.onUseRawPhaseImage)
    self.baselinePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.referencePhaseSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.tempMapSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.useThresholdFlagCheckBox.connect('toggled(bool)', self.onUseThreshold)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.baselinePhaseSelector.currentNode() and self.baselinePhaseSelector.currentNode() and self.tempMapSelector.currentNode()

  def onUseRawPhaseImage(self):
    pass
    
  def onUseThreshold(self):
    if self.useThresholdFlagCheckBox.checked == True:
      self.lowerThresholdSpinBox.enabled = True;      
      self.upperThresholdSpinBox.enabled = True;      
    else:
      self.lowerThresholdSpinBox.enabled = False;      
      self.upperThresholdSpinBox.enabled = False;      

  def onApplyButton(self):
    logic = PRFThermometryLogic()
    if self.useThresholdFlagCheckBox.checked == True:
      logic.run(self.rowPhaseImageFlagCheckBox.checked,
                self.baselinePhaseSelector.currentNode(), self.referencePhaseSelector.currentNode(),
                self.tempMapSelector.currentNode(), self.alphaSpinBox.value, self.gammaSpinBox.value,
                self.B0SpinBox.value, self.TESpinBox.value, self.BTSpinBox.value,
                self.upperThresholdSpinBox.value, self.lowerThresholdSpinBox.value)
    else:
      logic.run(self.rowPhaseImageFlagCheckBox.checked,
                self.baselinePhaseSelector.currentNode(), self.referencePhaseSelector.currentNode(),
                self.tempMapSelector.currentNode(), self.alphaSpinBox.value, self.gammaSpinBox.value,
                self.B0SpinBox.value, self.TESpinBox.value, self.BTSpinBox.value)

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

  def run(self, useRawPhaseImage, baselinePhaseVolumeNode, referencePhaseVolumeNode, tempMapVolumeNode, alpha, gamma, B0, TE, BT, upperThreshold=None, lowerThreshold=None):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(baselinePhaseVolumeNode, referencePhaseVolumeNode):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    print(baselinePhaseVolumeNode.GetName())
    print(referencePhaseVolumeNode.GetName())
    imageBaseline  = sitk.Cast(sitkUtils.PullFromSlicer(baselinePhaseVolumeNode.GetID()), sitk.sitkFloat64)
    imageReference = sitk.Cast(sitkUtils.PullFromSlicer(referencePhaseVolumeNode.GetID()), sitk.sitkFloat64)

    if tempMapVolumeNode:

      phaseDiff = None
      
      if useRawPhaseImage == True:
        imageBaselinePhase = imageBaseline*numpy.pi/4096.0
        imageBaselineReal = sitk.Cos(imageBaselinePhase)
        imageBaselineImag = sitk.Sin(imageBaselinePhase)
        imageBaselineComplex = sitk.RealAndImaginaryToComplex(imageBaselineReal, imageBaselineImag)
        
        imageReferencePhase = imageReference*numpy.pi/4096.0
        imageReferenceReal = sitk.Cos(imageReferencePhase)
        imageReferenceImag = sitk.Sin(imageReferencePhase)
        imageReferenceComplex = sitk.RealAndImaginaryToComplex(imageReferenceReal, imageReferenceImag)

        rotComplex = sitk.Divide(imageReferenceComplex, imageBaselineComplex)

        phaseDiff = sitk.ComplexToPhase(rotComplex)
        
      else:
        #imageTemp = (imageReference*2.0*numpy.pi/4096.0 - imageBaseline*2.0*numpy.pi/4096.0) / (alpha * gamma * B0 * TE) + BT
        print("(alpha, gamma, B0, TE, TE) = (%f, %f, %f, %f, %f)" % (alpha, gamma, B0, TE, BT))
        ## NOTE: gamma is given as Gyromagnetic ration / (2*PI) and needs to be mulitiplied by 2*PI
        phaseDiff = imageReference - imageBaseline
        
      imageTemp = (phaseDiff) / (alpha * 2.0 * numpy.pi * gamma * B0 * TE) + BT
        
      if upperThreshold or lowerThreshold:
        imageTempThreshold = sitk.Threshold(imageTemp, lowerThreshold, upperThreshold, 0.0)
        sitkUtils.PushToSlicer(imageTempThreshold, tempMapVolumeNode.GetName(), 0, True)
      else:
        sitkUtils.PushToSlicer(imageTemp, tempMapVolumeNode.GetName(), 0, True)

    logging.info('Processing completed')

    return True


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
