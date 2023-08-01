import numpy
from numpy.fft import fft2, fftshift
import matplotlib.pyplot as plot

from PIL import Image, ImageOps
from skimage.transform import (hough_line, hough_line_peaks)

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


#calculate assymmetry of hits accross detector

#@title
def ImageAsym( evtImage ):
  evtImage = evtImage.convert(mode="L")
  evtImage = ImageOps.invert(evtImage)
  theArray = numpy.asarray(evtImage, dtype = numpy.float32)
  h,w = numpy.shape(theArray)
  halfH = int(numpy.fix(0.5*h))
  halfW = int(numpy.fix(0.5*w))
  r = halfH-5
  box1 = theArray[halfH-r:halfH+r,0:2*r]
  box2 = theArray[halfH-r:halfH+r,halfW-r:halfW+r]
  box3 = theArray[halfH-r:halfH+r,w-2*r:w]
  frac1 = numpy.sum(box1) / numpy.sum(theArray)
  frac2 = numpy.sum(box2) / numpy.sum(theArray)
  frac3 = numpy.sum(box3) / numpy.sum(theArray)
  return numpy.std( (frac1, frac2, frac3) )

#Calculate occupancy
def ImageOccupancy( evtImage ):
  evtImage = evtImage.convert(mode="L")
  evtImage = ImageOps.invert(evtImage)
  theArray = numpy.asarray(evtImage, dtype = numpy.float32)
  total_sum = numpy.sum(theArray)
  max_sum = theArray.size*255
  percentage_fraction = total_sum/max_sum 
  #print("theAraray = ", theArray)
  #print(numpy.max(theArray))
  return percentage_fraction

#process XZ and YZ views for track angles

 #@title
def GetAnglesDet( evtImageXZ, evtImageYZ ):
  evtImageXZ = evtImageXZ.convert(mode="L")
  evtImageXZ = ImageOps.invert(evtImageXZ)
  evtImageYZ = evtImageYZ.convert(mode="L")
  evtImageYZ = ImageOps.invert(evtImageYZ)

  arrayXZ = numpy.asarray(evtImageXZ, dtype = numpy.float32)
  arrayYZ = numpy.asarray(evtImageYZ, dtype = numpy.float32)

  # Process XZ (top-down) view for azimuth angle
  hXZ,wXZ = numpy.shape(arrayXZ)
  halfH_XZ = int(numpy.fix(0.5*hXZ))
  halfW_XZ = int(numpy.fix(0.5*wXZ))

  r = halfH_XZ-5
  boxXZ = arrayXZ[halfH_XZ-r:halfH_XZ+r,halfW_XZ-r:halfW_XZ+r]

  fftXZ = fft2(boxXZ)/(wXZ*hXZ)
  fftXZ = fftshift(fftXZ,axes=(0,1))
  fftXZ = numpy.square(1000*numpy.abs(fftXZ)) #1000*numpy.abs(fftXZ)
  fftXZ = numpy.log(fftXZ)
  fftXZ[fftXZ < 8] = 0
  
  hspace, angles, distances = hough_line(fftXZ,theta=numpy.linspace(-numpy.pi/2,numpy.pi/2,3600))
  peakSpace, peakAngles, peakDists = hough_line_peaks(hspace, angles, distances, num_peaks = 10, min_angle = 1, min_distance = 1)
  peakAngles = peakAngles * (180./numpy.pi)
  goodAngles = []

  for angle in peakAngles:
    #if ( abs(  0 - abs(angle) ) < 1.0 ): continue
    if ( abs( 45 - abs(angle) ) < 1.0 ): continue
    if ( abs( 90 - abs(angle) ) < 1.0 ): continue
    if ( abs(135 - abs(angle) ) < 1.0 ): continue
    if ( abs(180 - abs(angle) ) < 1.0 ): continue
    if ( abs(225 - abs(angle) ) < 1.0 ): continue
    if ( abs(270 - abs(angle) ) < 1.0 ): continue
    if ( abs(315 - abs(angle) ) < 1.0 ): continue
    goodAngles.append(angle)
  
  total_angles_xz = len(goodAngles)
  if len(goodAngles) == 0: 
    return (0, 0, 0, 0)
  aveAngleXZ = numpy.average( goodAngles )
  stdAngleXZ = numpy.std( goodAngles )

  # Process YZ (side) view for calculating altitude angle.
  hYZ,wYZ = numpy.shape(arrayYZ)
  halfH_YZ = int(numpy.fix(0.5*hYZ))
  halfW_YZ = int(numpy.fix(0.5*wYZ))

  r = halfH_YZ-5
  boxYZ = arrayYZ[halfH_YZ-r:halfH_YZ+r,halfW_YZ-r:halfW_YZ+r]

  fftYZ = fft2(boxYZ)/(wYZ*hYZ)
  fftYZ = fftshift(fftYZ,axes=(0,1))
  fftYZ = numpy.square(1000*numpy.abs(fftYZ)) #1000*numpy.abs(fftYZ)
  fftYZ = numpy.log(fftYZ)
  fftYZ[fftYZ < 8] = 0

  hspace, angles, distances = hough_line(fftYZ,theta=numpy.linspace(-numpy.pi/2,numpy.pi/2,3600))
  peakSpace, peakAngles, peakDists = hough_line_peaks(hspace, angles, distances, num_peaks = 10, min_angle = 1, min_distance = 1)
  peakAngles = peakAngles * (180./numpy.pi)
  goodAngles = []
  
  for angle in peakAngles:
    #if ( abs(  0 - abs(angle) ) < 1.0 ): continue
    if ( abs( 45 - abs(angle) ) < 1.0 ): continue
    if ( abs( 90 - abs(angle) ) < 1.0 ): continue
    if ( abs(135 - abs(angle) ) < 1.0 ): continue
    if ( abs(180 - abs(angle) ) < 1.0 ): continue
    if ( abs(225 - abs(angle) ) < 1.0 ): continue
    if ( abs(270 - abs(angle) ) < 1.0 ): continue
    if ( abs(315 - abs(angle) ) < 1.0 ): continue
    goodAngles.append(angle)
  
  if len(goodAngles) == 0: 
    return (0, 0, 0, 0)
  aveAngleYZ = numpy.average( goodAngles )
  stdAngleYZ = numpy.std( goodAngles )

  # Find YZ angle from Hough lines on FFT
  trkAngleYZ = 90
  if aveAngleYZ < 0: trkAngleYZ = aveAngleYZ + 90
  else: trkAngleYZ = aveAngleYZ - 90

  # Find XZ angle from Hough lines on FFT and direction of tracks in YZ view
  trkAngleXZ = 0
  if trkAngleYZ > 0: trkAngleXZ = aveAngleXZ
  else: trkAngleXZ = 180 + aveAngleXZ

  if trkAngleXZ < 0: trkAngleXZ = trkAngleXZ + 360

  return trkAngleXZ, stdAngleXZ, trkAngleYZ, stdAngleYZ
  
#get altitude and azimuth of shower origin  

  #@title
def GetAltAz( trkAngleXZ, trkAngleYZ ):

  altAngle = numpy.abs( 180/numpy.pi * numpy.arctan( numpy.tan( numpy.pi/180. * (90. - trkAngleYZ)) * numpy.cos(trkAngleXZ*numpy.pi/180.) ) )

  # Chris Backhouse says docdb 5485 / numix 17 says FD is oriented -27o51'26'' from true North
  # or 332o03'58.071769" (clockwise from North). Email from Virgil Bocean to Alec Habig.
  FDrot = -27.857222 # degrees from North (along +z axis)
  azAngle = trkAngleXZ + FDrot
  if azAngle < 0: azAngle = azAngle + 360

  return altAngle, azAngle
  
#process NOvA event timestamp

  #@title
def ProcessTimestamp( evtTimestamp ):
  evtNanosecs = evtTimestamp & 0xFFFFFFFF
  evtTimestamp = (evtTimestamp >> 32) & 0xFFFFFFFF
  evtTime = Time(evtTimestamp + evtNanosecs*1e-9,format='unix') 
  return evtTime

#create sky coordinate object
  
  #@title
def MakeSkyCoord( altAngle, azAngle, evtTime ):

  # FD position in degrees.
  # These come from clicking on the FD building in Google Maps.
  FDlat  =  48.378574
  FDlon = -92.831051
  # Another number we've used is 354.9. PM with Alex. Weather station?
  FDheight = 367.6 # meters. docdb 1206 p8

  FDloc = EarthLocation(lat=FDlat*u.deg, lon=FDlon*u.deg, height=FDheight*u.m)

  theAltAzCoord = AltAz(alt = altAngle*u.deg, az = azAngle*u.deg, obstime = evtTime, location = FDloc)

  return SkyCoord(theAltAzCoord)

  
  

