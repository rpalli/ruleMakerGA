#!/usr/bin/env python

import numpy, ctypes

class MyImage:
  def __init__(self):
    self.lib = ctypes.cdll.LoadLibrary('./testRun.so')
    self.nrows = self.lib.image_nrows()
    self.ncols = self.lib.image_ncols()
    self.data = numpy.empty((self.nrows,self.ncols), dtype=numpy.double)
    self.data2= numpy.zeros((self.nrows,self.ncols), dtype=numpy.int8)
  def get(self):
    self.lib.image_get(ctypes.c_void_p(self.data.ctypes.data),ctypes.c_void_p(self.data2.ctypes.data) )

if __name__ == '__main__':
  updateBooler=ctypes.cdll.LoadLibrary('./testRun.so')
  updateBool=updateBooler.updateBool  
  #img = MyImage()
  #print(img.nrows,img.ncols)
  #print(img.data)
  #img.get()
  currentNode=0
  oldValue = numpy.zeros(2, dtype=numpy.int8)
  oldValue[1]=1
  print(currentNode)
  nodeIndividual=numpy.zeros(2, dtype=numpy.int8)
  print(nodeIndividual)
  indLen=1
  andLenList= numpy.zeros(300, dtype=numpy.int8)
  andLenList[0]=1
  andNodes=numpy.zeros([300,128],  dtype=numpy.int8)
  andNodes[0][0]=1
  andNodeInvertList=numpy.array([[False for i in range(128)] for j in range(300)]) 
  
  currentNode=ctypes.c_void_p(currentNode)

  nodeIndividual=ctypes.c_void_p(nodeIndividual.ctypes.data)
  indLen=ctypes.c_void_p(indLen)
  andNodes=ctypes.c_void_p(andNodes.ctypes.data)
  andLenList=ctypes.c_void_p(andLenList.ctypes.data)
  oldValue=ctypes.c_void_p(oldValue.ctypes.data)
  andNodeInvertList=ctypes.c_void_p(andNodeInvertList.ctypes.data)
  value=updateBool(currentNode,oldValue,nodeIndividual, indLen,andNodes, andLenList,andNodeInvertList)
  print(value)
  #print(img.data)
