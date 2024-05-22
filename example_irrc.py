import sys
import os
import glob
from astropy.io import fits
import numpy as np
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe

def irrc_correct_image_cube(dataIn, calFile,multiThread:bool=True, externalPixelFlags:np.ndarray=None, superbias_corrected:bool=True):
    """Returns image cube corrected using IRRC in memory, given calibration file."""

    #Grab coefficients from calFile
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)
    #De-stripe input
    dataIn_tmp = destripe(dataIn)
    #Apply to image cube
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags)
    return dataIRRC

def irrc_correct_frame(dataIn, superbias, calFile,multiThread:bool=True, externalPixelFlags:np.ndarray=None, superbias_corrected:bool=True):
    """Returns image cube corrected using IRRC in memory, given calibration file."""

    #Grab coefficients from calFile
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)
    #De-stripe input (subtract off superbias)
    dataIn_tmp = destripe(dataIn[None,:,:])
    dataIn_tmp[:,:,:4096] = dataIn_tmp[:,:,:4096]-superbias[None,:,:]
    #Apply to image cube (after correction add superbias back in)
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags, superbias_corrected=True)+superbias
    return dataIRRC

if __name__ == "__main__":
    image_hdu = 1
    sca = 'C1'

    if sca not in ['C1','C2','C3','C4']:
        print("Please use chip name C1, C2, C3, or C4. Exiting.")
        sys.exit(1)

    # libdir = r'C:\PychcarmProjects\PRIMErampfit\example_irrc'
    libdir = os.path.join(os.path.dirname(__file__), 'IRRC_calfiles')
    calFile = os.path.join(libdir,"irrc_weights_"+sca+".h5")

    #Provide bad reference pixel mask
    maskFile = os.path.join(libdir, sca+"_bad_ref_pix_mask.fits")
    print("Reading bad reference pixel mask: %s"%(maskFile))
    mask = fits.open(maskFile)[1].data>0

    #Provide superbias file
    superbiasFile = glob.glob(os.path.join(libdir, 'super_bias'+sca+'*'))[0]
    print("Using superbias: %s"%(superbiasFile))
    superbias = fits.open(superbiasFile)[0].data
    #Mask bad reference pixels
    supercpy = superbias.copy()
    supercpy[mask[:,:4096]] = 0

    #Read frame list
    data_dir = r'D:\NLC\{}'.format(sca)
    frame_list = ['01124975C1.fits.fz']
    frame_list = [os.path.join(data_dir, f) for f in frame_list]

    #Form image cube in memory
    img_cube_list = []
    hdr_list = []
    for f in frame_list:
        hdu =  fits.open(f)
        d = hdu[image_hdu].data[:,6:] #Assumes data is in 1st HDU of FITS and headers are present
        h = hdu[image_hdu].header
        hdr_list.append(h)
        #Mask bad reference pixels
        dcpy = d.copy()
        dcpy[mask] = 0
        img_cube_list.append(dcpy)
    dataIn = np.asarray(img_cube_list)

    #Flag to correct using cube in memory
    correct_cube = False
    if correct_cube:
        #Apply IRRC correction
        print("Applying IRRC correction to image cube in memory")
        dataOut = irrc_correct_image_cube(dataIn, calFile, superbias_corrected=False)

        #Write output
        for i,f in enumerate(frame_list):
            fout = f.split('.fits')[0]+'_ircc.fits'
            print("Writing output to %s"%(fout))
            outhdu = fits.PrimaryHDU(dataOut[i],header=hdr_list[i])
            outhdul = fits.HDUList([outhdu])
            outhdul.writeto(fout, overwrite=True)

    #Flag to correct by frame
    correct_frame = True
    if correct_frame:
        #Apply IRRC correction
        print("Applying IRRC correction to image frames in memory")
        #Write output
        for i,f in enumerate(frame_list):
            fout = f.split('.fits')[0]+'_ircc.fits'
            print("Writing output to %s"%(fout))
            dataOut = irrc_correct_frame(dataIn[i], supercpy, calFile, superbias_corrected=False)
            outhdu = fits.PrimaryHDU(dataOut,header=hdr_list[i])
            outhdul = fits.HDUList([outhdu])
            outhdul.writeto(fout, overwrite=True)


