#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys,os
import string
import math

import argparse

def read_plugmap(filename):
    debug=False
    file=open(filename,"r")
    doc={}
    intypedef=False
    indices={}
    indices["HOLETYPE"]=8
    indices["OBJECT"]=21
    indices["ra"]=9
    indices["dec"]=10
    indices["xfoc"]=22
    indices["yfoc"]=23
    objects={}
    for k in indices :
        objects[k]=[]
        
    for line in file.readlines() :
        line=line.strip().replace('\t',' ')
        if debug :
            print "line: ",line
        if len(line)==0 :
            continue
        if line[0]=="#" :
            continue
        if line.find("typedef")>=0 :
            intypedef=True
            if debug :
                print "now in typedef"
            continue
        if intypedef and line.find("}")>=0 :
            intypedef=False
            if debug :
                print "end of typedef"
            continue
        if intypedef :
            continue
        
        if line.find("PLUGMAPOBJ")>=0 :
            
            tmp=line.split(" ")
            entries=[]
            for t in tmp :
                if len(t)>0 :
                    entries.append(t)
            
            for k in objects.keys() :
                i=indices[k]
                val=entries[i]
                #print k,i,val
                tmp=None
                try : 
                    tmp=string.atoi(val)
                except ValueError :
                    pass
                if tmp is None :
                    try : 
                        val=string.atof(val)
                    except ValueError :
                        pass
                if tmp is not None :
                    val=tmp
                objects[k].append(val)
                if debug :
                    print "added one PLUGMAPOBJ" 
            continue

        tmp=line.strip().split(" ")
        entries=[]
        for t in tmp :
            if len(t)>0 :
                entries.append(t)
        
        if len(entries)>=2 :
            key=entries[0]
            val=entries[1]
            tmp=None
            try : 
               tmp=string.atoi(val)
            except ValueError :
                pass
            if tmp is None :
                try : 
                    val=string.atof(val)
                except ValueError :
                    pass
            if tmp is not None :
                val=tmp
            doc[key]=val
            if debug :
                print "added doc",key,val 

    # convert objects into np.array
    for k in objects :
        objects[k]=np.array(objects[k])
    return doc,objects


class OpticalDistortion() :
    def __init__(self,platescale) :

        self.platescale=platescale # has units
        
        # see ~/software/platedesign/trunk/pro/plate/ad2xyfocal.pro
        coef=np.array([-0.000137627, -0.00125238, 1.5447e-09, 
                        8.23673e-08, -2.74584e-13, -1.53239e-12, 
                        6.04194e-18, 1.38033e-17, -2.97064e-23, 
                        -3.58767e-23])
        self.achromatic_distortion_pol=np.poly1d(coef[::-1])

        # see ~/software/platedesign/trunk/pro/plate/apo_rdistort.pro
        mm_per_rad =platescale*180/math.pi
        self.chromatic_distort_radii=np.arcsin(np.linspace(0,90,10)*math.pi/(60*180))*mm_per_rad
        print "RADII=",self.chromatic_distort_radii
        
        self.chromatic_distort_wave=np.array([5300,4000,5500,6000,8000,10000,15350,15950,16550])
        nw=self.chromatic_distort_wave.size
        nr=self.chromatic_distort_radii.size
        
        self.chromatic_distort=np.array([
                [0.,36.26,72.53,108.84,145.18,181.53,217.90,254.29,290.77,327.44],
                [0.,-0.002,-0.003,-0.004,-0.005,-0.005,-0.005,-0.004,-0.002,0.003],
                [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                [0.,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,-0.001],
                [0.,0.001,0.003,0.003,0.004,0.004,0.004,0.003,0.002,-0.003],
                [0.,0.002,0.004,0.005,0.005,0.005,0.005,0.005,0.003,-0.004],
                [0.,0.003,0.006,0.007,0.008,0.008,0.008,0.008,0.004,-0.006],
                [0.,0.003,0.006,0.008,0.008,0.009,0.009,0.008,0.004,-0.006],
                [0.,0.004,0.006,0.008,0.009,0.009,0.009,0.008,0.004,-0.007]])
        
        # apply scaling
        scale=np.zeros((nr))
        scale[1:]=self.chromatic_distort_radii[1:]/self.chromatic_distort[0,1:]
        self.chromatic_distort[1:] *= scale
        self.chromatic_distort[0]=0.
        # sort wave
        ii=np.argsort(self.chromatic_distort_wave)
        self.chromatic_distort_wave=self.chromatic_distort_wave[ii]
        for j in range(nr) :
            self.chromatic_distort[:,j]=self.chromatic_distort[ii,j]
        # in ad2xyfocal, a reference wavelength of 5000A instead of 5500A is used !!
        ref_distort = np.zeros((nr))
        for j in range(nr) :
            ref_distort[j]=np.interp(5000,self.chromatic_distort_wave,self.chromatic_distort[:,j])
        self.chromatic_distort -= ref_distort
        
        """
        plt.plot(self.chromatic_distort_wave,self.chromatic_distort[:,-1],"o-")
        
        ww=np.linspace(4000,8000,200)*u.angstrom
        r=self.chromatic_distort_radii[-1]
        dd=np.zeros((ww.size))
        for i in range(ww.size) :
            dd[i]=self.chromatic_distortion(r,ww[i]).to(u.mm).value
        plt.plot(ww,dd,c="r")
        plt.show()
        """


    def chromatic_distortion(self,radius,wavelength) : # with radius and wave with units , returns delta r to be added
        i=np.where(self.chromatic_distort_wave>=wavelength)[0]
        if i.size == 0 :
            i=1
        else :
            i=min(max(1,i[0]),self.chromatic_distort_radii.size-1)
        dist1=np.interp(radius,self.chromatic_distort_radii,self.chromatic_distort[i-1])
        dist2=np.interp(radius,self.chromatic_distort_radii,self.chromatic_distort[i])
        dist=np.interp(wavelength,[self.chromatic_distort_wave[i-1],self.chromatic_distort_wave[i]],[dist1,dist2])
        return dist
    
    def distortion(self,radius,wavelength) : 
        return self.achromatic_distortion_pol(radius) + self.chromatic_distortion(radius,wavelength)


# same result as idlutils/goddard/pro/astro/hadec2altaz.pro 
# but with adr calibrated using astropy
def hadec2altaz(ha, dec, lat, wavelength=None) : # ha,dec,lat in deg, wave in a, returns alt,az
    d2r = math.pi/180.
    sh = math.sin(ha*d2r)
    ch = math.cos(ha*d2r)
    sd = math.sin(dec*d2r)
    cd = math.cos(dec*d2r)
    sl = math.sin(lat*d2r)
    cl = math.cos(lat*d2r)
    
    """
    x=np.array([cd*ch,cd*sh,sd])
    r=np.array([[sl,0,-cl],[0,1,0],[cl,0,sl]])
    x=r.dot(x)    
    x0=x[0]
    x1=x[1]
    x2=x[2]
    """
    x0 = - ch * cd * sl + sd * cl
    x1 = - sh * cd
    x2 = ch * cd * cl + sd * sl
    
    r=math.sqrt(x0**2+x1**2)
    az = math.atan2(-x1,-x0) /d2r
    alt = math.atan2(x2,r) / d2r
    
    if wavelength is not None :
        # arcsec per unit of tan(zenith)
        fact=np.interp(wavelength,[3000,3500,4000,5000,5400,6000,7000,8000],[44.166347,43.365612,42.8640697818,42.292551282,42.1507465805,41.990386,41.811009,41.695723])
        alt += fact*(r/x2)/3600.
    
    return alt,az

# exact same routine as altaz2rpa in idl, needed to get same platescale definition
def altaz2xy(alt,az,altcen,azcen,platescale) :
    d2r=math.pi/180
    xx= -np.sin(az*d2r) * np.sin((90-alt)*d2r)
    yy= -np.cos(az*d2r) * np.sin((90-alt)*d2r)
    zz= np.cos((90-alt)*d2r)
    xi= -xx*np.cos(azcen*d2r) + yy*np.sin(azcen*d2r)
    yi= -yy*np.cos(azcen*d2r) - xx*np.sin(azcen*d2r)
    zi= zz
    xl= xi
    yl= yi*np.sin((90-altcen)*d2r) + zi*np.cos((90-altcen)*d2r)
    zl= zi*np.sin((90-altcen)*d2r) - yi*np.cos((90-altcen)*d2r)
    
    rfocal=np.arcsin(np.sqrt(xl**2+zl**2))/d2r*platescale
    posang=np.arctan2(-xl, zl)
    return rfocal*np.cos(posang),rfocal*np.sin(posang)


def hadec2xy(ha,dec,alt0,az0,crot,srot,latitude,platescale,distortion,wavelength) :
    alt,az = hadec2altaz(ha,dec,latitude,wavelength)
    x,y    = altaz2xy(alt,az,alt0,az0,platescale)
    rscale = 1
    if 1 :
        # Distortion, see ad2xyfocal.pro
        r = np.sqrt(x**2 + y**2)
        if r>0 :
            rscale = 1+distortion.distortion(r,wavelength)/r
    
    
    # Rotate the focal plane so that +y points towards a point that is offset from
    # the plate center along DEC by +1.5 degrees.        
    xr = rscale*(x*crot-y*srot)
    yr = rscale*(x*srot+y*crot)
    return -xr,yr,alt,az
        
def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', type=str, default='plPlugMapP-4392.par',
        help='Input plugmap filename.')
    parser.add_argument('--output', type=str, default='myPlugMapP-4392.list',
        help='Output filename.')
    parser.add_argument('--ha', type=float, default=0,
        help='Design hour angle (degrees).')
    args = parser.parse_args()

    filename = args.input
    ofilename = args.output
    ha_obs = args.ha

    doc, objects = read_plugmap(filename)

    ra=objects["ra"]
    dec=objects["dec"]
    xfoc=objects["xfoc"]
    yfoc=objects["yfoc"]

    ha_design=doc["ha"]
    ra0=doc["raCen"]
    dec0=doc["decCen"]
    mjd=doc["mjdDesign"]
    print "design MJD=%d HA=%f ra=%f dec=%f"%(mjd,ha_design,ra0,dec0)

    # APO lat=32.7797556 in plate_refrac.pro
    latitude=32.7797556

    # optical distortion
    # from platedesign/trunk/pro/plate/get_platescale.pro
    platescale = 217.7358 
    distortion = OpticalDistortion(platescale)

    # only reference for this wavelength I could find is in code platedesign/trunk/pro/plate/adr.pro
    refwave=5400.0
    gal=np.where(objects["OBJECT"]=="GALAXY")[0]
    qso=np.where(objects["OBJECT"]=="QSO")[0]
    star=np.where(objects["OBJECT"]=="SPECTROPHOTO_STD")[0]
    na=np.where(objects["OBJECT"]=="NA")[0]
    nobj=xfoc.size
    wave_design=refwave*np.ones((nobj))
    wave_design[gal]=5400.
    wave_design[qso]=4000.
    wave_design[star]=5400.
    wave_obs=7450*np.ones((nobj))
    wave_obs[gal]=7450. # to study r1/r2
    wave_obs[qso]=7450.
    wave_obs[star]=7450.

    # for design
    alt0_design,az0_design     = hadec2altaz(ha_design, dec0, latitude, refwave)
    print "Design ALT (ref wave)=",alt0_design
    print "Design AZ  (ref wave)=",az0_design
    # rotation of plate to get vertical dec
    altfid,azfid = hadec2altaz(ha_design, dec0+1.5, latitude, refwave)
    xfid,yfid      = altaz2xy(altfid,azfid,alt0_design,az0_design,platescale)
    rotation_angle = np.arctan2(xfid,yfid)
    crot_design           = np.cos(rotation_angle)
    srot_design           = np.sin(rotation_angle)
    # same for obs
    alt0_obs,az0_obs     = hadec2altaz(ha_obs, dec0, latitude, refwave)
    print "Obs ALT (ref wave)=",alt0_obs
    print "Obs AZ  (ref wave)=",az0_obs
    # rotation of plate to get vertical dec
    altfid,azfid = hadec2altaz(ha_obs, dec0+1.5, latitude, refwave)
    xfid,yfid      = altaz2xy(altfid,azfid,alt0_obs,az0_obs,platescale)
    rotation_angle = np.arctan2(xfid,yfid)
    crot_obs           = np.cos(rotation_angle)
    srot_obs           = np.sin(rotation_angle)

    # compute, at design hour angle = ha_design
    xdesign=np.zeros((nobj))
    ydesign=np.zeros((nobj))
    alt_design=np.zeros((nobj))
    az_design=np.zeros((nobj))
    # compute, at observed hour angle = ha_obs
    xobs=np.zeros((nobj))
    yobs=np.zeros((nobj))
    alt_obs=np.zeros((nobj))
    az_obs=np.zeros((nobj))

    selection=range(nobj)

    for o in selection :
        x,y,alt,az = hadec2xy(ha_design-(ra[o]-ra0),dec[o],alt0_design,az0_design,crot_design,srot_design,latitude,platescale,distortion,wave_design[o])
        xdesign[o] = x
        ydesign[o] = y
        alt_design[o] = alt
        az_design[o] = az

        x,y,alt,az = hadec2xy(ha_obs-(ra[o]-ra0),dec[o],alt0_obs,az0_obs,crot_obs,srot_obs,latitude,platescale,distortion,wave_obs[o])
        xobs[o] = x
        yobs[o] = y
        alt_obs[o] = alt
        az_obs[o] = az

    file=open(ofilename,"w")
    file.write("#ra dec xfoc yfoc wavedesign xdesign ydesign altdesign azdesign waveobs xobs yobs altobs azobs hole obj\n")
    for o in selection :
        file.write("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s\n"%(ra[o],dec[o],xfoc[o],yfoc[o],wave_design[o],xdesign[o],ydesign[o],alt_design[o],az_design[o],wave_obs[o],xobs[o],yobs[o],alt_obs[o],az_obs[o],objects["HOLETYPE"][o],objects["OBJECT"][o]))
    file.close()
    print "wrote", ofilename

    
if __name__ == '__main__':
    main()

