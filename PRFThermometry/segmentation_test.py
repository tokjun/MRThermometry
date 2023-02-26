#!/usr/bin/env python3

import sitkUtils
import SimpleITK as sitk

def segment_test(baselineNode, referenceNode, outputName):

    baselineImage = sitk.Cast(sitkUtils.PullVolumeFromSlicer(baselineNode), sitk.sitkInt16)
    referenceImage = sitk.Cast(sitkUtils.PullVolumeFromSlicer(referenceNode), sitk.sitkInt16)

    # Subtraction
    image_diff = baselineImage - referenceImage

    # Threshold
    otsu_filter = sitk.MaximumEntropyThresholdImageFilter()
    otsu_filter.SetInsideValue(0)
    otsu_filter.SetOutsideValue(1)
    otsu_filter.SetNumberOfHistogramBins(5)
    image_threshold = otsu_filter.Execute(image_diff)

    # Relabel connected regions
    cc_filter = sitk.ConnectedComponentImageFilter()
    image_cc = cc_filter.Execute(image_threshold)
    relabel_filter = sitk.RelabelComponentImageFilter()
    relabel_filter.SetMinimumObjectSize(30)
    image_relabel = relabel_filter.Execute(image_cc)

    # Pick the largest region (regions are sorted by area)
    image_obj = (image_relabel == 1)

    # Smoothing
    dilate_filter = sitk.BinaryDilateImageFilter()
    dilate_filter.SetKernelRadius([1,1,1])
    dilate_filter.SetForegroundValue(1.0)
    dilate_filter.SetBackgroundValue(0.0)
    image_obj = dilate_filter.Execute(image_obj)

    erode_filter = sitk.BinaryErodeImageFilter()
    erode_filter.SetKernelRadius([2,2,2])
    erode_filter.SetForegroundValue(1.0)
    erode_filter.SetBackgroundValue(0.0)
    image_obj = erode_filter.Execute(image_obj)

    outputNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode')
    outputNode.SetName(outputName)
    sitkUtils.PushVolumeToSlicer(image_obj, outputName, 0, True)
