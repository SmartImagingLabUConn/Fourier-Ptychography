Content
-- Data
   -- blood_aberration.mat
   -- bloodsmear_blue.mat
   -- bloodsmear_green.mat
   -- bloodsmear_red.mat
   -- HE_blue.mat
   -- HE_green.mat
   -- HE_red.mat
   -- IHC_blue.mat
   -- IHC_green.mat
   -- IHC_red.mat
   -- MouseKidney_green.mat
   -- USAF_red.mat
All samples are captured using a Nikon 2X objective with NA = 0.1.
The pixel size of low-resolution raw images on the sample plane is 1.845 microns.
The slide glass thickness of each sample is 1 mm with refraction index RI = 1.52.

Illumination strategy:
A 32 by 32 LED array with a pitch of 4 mm is used in the setup.
The coordinate is for the first LED is (x=18,y=20).
15*15 LEDs are lit sequentially in a spiral-out manner.
The distance between the LED array and the sample is 90.88 mm.

Each MAT-file contains:
aberration: pre-calibrated aberration, if available
imlow_HDR: low-resolution measurements
theta: rotation angle of LED array to the camera sensor frame, in degree
wlength: central peak wavelength of LED, in m
xint,yint: offset of initial LED to the patch center, in mm
z: known defocus distance, in m
