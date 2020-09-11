"""
    Super Class of Geodetic Time Series Class & methods (Sgts)
    Sgts is a record of Gts and enables to apply methods at the same time to various Gts
        
"""

import sys
from pyacs.gts.Gts import Gts

class Sgts:


    ###################################################################
    def __init__ (self, ts_dir='.', add_key='', verbose=True, name_filter='', read=True,sites=[],lexclude=[],type=None,xyz=True):
    ###################################################################
        self.dir=ts_dir
        self.filter=name_filter
        self.verbose=verbose
        if read:
            self.read_ts(ts_dir=self.dir, name_filter=self.filter, add_key=add_key, verbose=verbose,sites=sites,lexclude=lexclude,type=type,xyz=xyz)



###################################################################
##  METHODS IMPORT
###################################################################

import pyacs.gts.Sgts_methods.add_offsets_dates
import pyacs.gts.Sgts_methods.append
import pyacs.gts.Sgts_methods.copy
import pyacs.gts.Sgts_methods.delts
import pyacs.gts.Sgts_methods.frame
import pyacs.gts.Sgts_methods.gts
import pyacs.gts.Sgts_methods.has_ts
import pyacs.gts.Sgts_methods.lGts
import pyacs.gts.Sgts_methods.n
import pyacs.gts.Sgts_methods.lcode
import pyacs.gts.Sgts_methods.medvel
import pyacs.gts.Sgts_methods.read_gmt
import pyacs.gts.Sgts_methods.read_gts_conf
import pyacs.gts.Sgts_methods.read_soln
import pyacs.gts.Sgts_methods.read_ts
import pyacs.gts.Sgts_methods.same_site
import pyacs.gts.Sgts_methods.save_velocity
import pyacs.gts.Sgts_methods.sel_rectangle
import pyacs.gts.Sgts_methods.sel_period
import pyacs.gts.Sgts_methods.show_map
import pyacs.gts.Sgts_methods.stat_site
import pyacs.gts.Sgts_methods.sub
import pyacs.gts.Sgts_methods.to_displacement
import pyacs.gts.Sgts_methods.write_pck
import pyacs.gts.Sgts_methods.common_mode


Sgts.append             = pyacs.gts.Sgts_methods.append.append
Sgts.copy               = pyacs.gts.Sgts_methods.copy.copy
Sgts.delts              = pyacs.gts.Sgts_methods.delts.delts
Sgts.frame              = pyacs.gts.Sgts_methods.frame.frame
Sgts.gts                = pyacs.gts.Sgts_methods.gts.gts
Sgts.has_ts             = pyacs.gts.Sgts_methods.has_ts.has_ts
Sgts.lGts               = pyacs.gts.Sgts_methods.lGts.lGts
Sgts.n                  = pyacs.gts.Sgts_methods.n.n
Sgts.lcode              = pyacs.gts.Sgts_methods.lcode.lcode
Sgts.medvel             = pyacs.gts.Sgts_methods.medvel.medvel
Sgts.read_gmt           = pyacs.gts.Sgts_methods.read_gmt.read_gmt
Sgts.read_gts_conf      = pyacs.gts.Sgts_methods.read_gts_conf.read_gts_conf
Sgts.read_soln          = pyacs.gts.Sgts_methods.read_soln.read_soln
Sgts.read_ts            = pyacs.gts.Sgts_methods.read_ts.read_ts
Sgts.same_site          = pyacs.gts.Sgts_methods.same_site.same_site
Sgts.save_velocity      = pyacs.gts.Sgts_methods.save_velocity.save_velocity
Sgts.sel_rectangle      = pyacs.gts.Sgts_methods.sel_rectangle.sel_rectangle
Sgts.sel_period         = pyacs.gts.Sgts_methods.sel_period.sel_period
Sgts.show_map           = pyacs.gts.Sgts_methods.show_map.show_map
Sgts.stat_site          = pyacs.gts.Sgts_methods.stat_site.stat_site
Sgts.sub                = pyacs.gts.Sgts_methods.sub.sub
Sgts.to_displacement    = pyacs.gts.Sgts_methods.to_displacement.to_displacement
Sgts.write_pck          = pyacs.gts.Sgts_methods.write_pck.write_pck
Sgts.common_mode        = pyacs.gts.Sgts_methods.common_mode.common_mode


# 
#     ###################################################################
#     def gts(self , method , *args , **kwargs):
#     ###################################################################
#         """
#         apply a gts method to all Gts instance of the current Sgts object
#          
#         :param method: Gts method to be applied as string
#         :param *arg: arguments for the Gts method to be applied 
#         :param **kwarg: keyword arguments for the Gts method to be applied 
#         
#         :example : ts.gts('detrend',periods=[2010.0,2013.0])
#         """
# 
#         verbose = kwargs.get('verbose', False)        
#          
#         new_ts = Sgts(read=False)
#          
#         lsite=self.lcode()
#          
#         for site in sorted(lsite):
#             if verbose:
#                 print('-- processing ', site)
#             
#             try:
#                 func = getattr(self.__dict__[site], method)
#                 new_ts.append(func(*args, **kwargs) )
#             except:
#                 print('!!!WARNING: problem with method %s on gts %s. Removed from output' % ( method, site ))
#         return( new_ts )
# 
# 
#     ###################################################################
#     def stack(self , verbose=True):
#     ###################################################################
#         """
#         Compute the stack of the Gts in the current Sgts
#          
#         :return Gts: the stacked Gts with '_STK' code
#         
#         :note 1: the stack is performed on the .data attribute only with uncertainties set to zero.
#         :note 2: No detrend is applied.
#         :note 3: Gap in data are not handled and time series are assumed to have exactly the same number of data at the same dates.
#         """
#         
#         
#         
#         new_gts = Gts()
#         new_gts.code = '_STK'
#         
#         data = None
#         
#         for code in self.lcode():
#             if verbose:
#                 print("-- adding %s in stack." % code)
#             
#             wgts = self.__dict__[code]
#             try:
#                 if data is None:
#                     data = wgts.data
#                     np_date = wgts.data[:,0]
#                 else:
#                     data = data + wgts.data
#             except:
#                 print("!!!WARNING: some problem (likely missing data) for time series: %s" % code)
#             
#         new_gts.data = data
#         new_gts.data[:,0] = np_date
#                 
#         return new_gts
# 
#         
# 
# ###################################################################
#     def lcode(self,lexclude=[],linclude=[]):
# ###################################################################
#         """
#         Returns the list of Gts codes in the current Sgts
#         exclude is a list of code to be excluded
#         
#         :param lexclude: list of sites to be excluded
#         :param linclude: list of sites to be included, excluding all other.
#         
#         """
#         # import
#         import re 
#         
#         
#         lcode=[]
#         if linclude != []:
#             for k in list(self.__dict__.keys()):
#                 if ((re.match('[A-Z0-9_]{4}',k))) and (k not in lexclude) and (k in linclude):lcode.append(k)
#         else:
#             for k in list(self.__dict__.keys()):
#                 if ((re.match('[A-Z0-9_]{4}',k))) and (k not in lexclude):lcode.append(k)
# 
#         return(lcode)
# 
# ###################################################################
#     def lGts(self,lexclude=[],linclude=[]):
# ###################################################################
#         """
#         Returns the list of Gts in the current Sgts
#         :param lexclude: list of sites to be excluded
#         :param linclude: list of sites to be included, excluding all other.
#         
#         """
#         
#         lsite=self.lcode(lexclude=lexclude,linclude=linclude)
#         Lgts=[]
#         for site in sorted(lsite):Lgts.append(self.__dict__[site])
#         return(Lgts)
# 
# 
# ###################################################################
#     def delts(self,code):
# ###################################################################
#         """
#         Delete a time series from an Sgts instance
# 
#         :param code: code to be excluded
#         """
# 
#         delattr(self, code)
# 
# 
# ###################################################################
#     def sub(self,lexclude=[],linclude=[]):
# ###################################################################
#         """
#         Returns a new Sgts instance excluding Gts with code in lexclude and keeping Gts with code in include
# 
#         :param lexclude: list of sites to be excluded
#         :param linclude: list of sites to be included, excluding all other.
#         
#         """
#         
#             
#         if linclude != []:
#             sub_Sgts = Sgts(read=False)
#             for code in linclude:
#                 if code not in lexclude:
#                     sub_Sgts.append(self.__dict__[code].copy())
#         else:
#             sub_Sgts = self.copy()
#             for code in lexclude:
#                 sub_Sgts.delts(code)
#         
#         return(sub_Sgts)
# 
# 
# ###################################################################
#     def append(self, gts):
# ###################################################################
#         """
#         Append a Gts to the current Sgts
#         
#         :param gts: Gts instance to be appended to the current Sgts instance
#         """
#         
#         self.__dict__[gts.code]=gts
# 
# ###################################################################
#     def has_ts(self, code):
# ###################################################################
#         """
#         Tests whether a time series exists in the current Sgts instance
#         
#         :param code: 4-character code
#         """
#         
#         def has_key(key, data):
#             return True if key in data else False
# 
#         if has_key(code , self.__dict__):
#             return True
#         else:
#             return False
# 
# 
# ###################################################################
#     def common_mode(self, lGts, lGts_ref, date, verbose=True):
# ###################################################################
#         """
#         Common mode filter
#         
#         :todo use Sgts objects instead of lists
#         
#         lGts    : list of time series to which the common mode will be applied
#         lGts_ref: list of time series used to calculates the common mode
#         date    : period where common mode will be applied
#         Output:
#         lGts    : filtered time series
#         common_ts: time series of common mode
#         """
#         
#         # import
#         import numpy as np
#         import pyacs.lib.astrotime
# 
#         # homogenize dates
#         for gts in lGts_ref:
#             np_dates=gts.data[:,0]
#             np_dates_mjd = pyacs.lib.astrotime.decyear2mjd(np_dates)
#             np_dates_mjd = np.trunc( np_dates_mjd ) + 0.5
#             gts.data[:,0] = np_dates_mjd
# 
#         # converts mjd dates back to decimal year
#         for gts in lGts_ref:
#             gts.data[:,0] = pyacs.lib.astrotime.mjd2decyear(gts.data[:,0])
# 
# 
#         # get the list of dates of ref
#         ldates=np.array([])
#         for gts in lGts_ref:
#             lindex=np.where((gts.data[:,0] >= date[0]) & (gts.data[:,0] <= date[-1]))
#             ldates=np.append(ldates,gts.data[lindex,0])
#             ldates=np.unique(np.array(ldates))
#         # calculates the common mode
#         if verbose:print("-- Calculating common mode")
#         common_mode=Gts(code='CMM1')
#         common_mode.data=np.zeros((1,7))
# 
#         for date in ldates:
# 
#             common_mode_date=np.zeros((1,7))
#             for gts in lGts_ref:
#                 #print "  -- processing ",gts.code
#                 if date in list(gts.data[:,0]):
#                     index=np.where(date==gts.data[:,0])
#                     common_mode_date=np.vstack((common_mode_date,gts.data[index[0],:7]))
#                 else:
#                     print('-- date missing for ',gts.code)
#             common_mode_date=np.delete(common_mode_date,0,axis=0)
#             
#             common_mode.data=np.vstack((common_mode.data,np.mean(common_mode_date,axis=0)))
#         common_mode.data=np.delete(common_mode.data,0,axis=0)
#         
#         common_mode.write_mb_file('.')  
# 
#         if verbose:print("-- Removing common mode")
#         
#         lfiltered_gts=[]
#         
#         for gts in lGts:
#             print("  -- ",gts.code)
#             filtered_ts=gts.substract_ts(common_mode)
#             if filtered_ts != None:
#                 lfiltered_gts.append(filtered_ts)
#             else:
#                 print(gts.code," has no data common with reference sites")
#         
#         
# #        common_mode.data[:,0] = pyacs.mjd2decyear(common_mode.data[:,0])
#         
#         
#         # returns results
#         return(lfiltered_gts,common_mode)
# 
# ###################################################################
#     def write_pck(self,outfile, verbose=True):
# ###################################################################
#         """
#         writes a Sgts object as a pck (pickle)
#         
#         :param outfile: output file name. If not provided, a pck extension will be added.
#         :param verbose: verbose mode
#         
#         """
#         # import
#         import pickle
#         import os
#     
#         # add pck extension if not provided
#         
#         if outfile[-4:] != '.pck':
#             outfile = outfile+'.pck'
# 
#         os.makedirs( os.path.dirname( outfile ) , exist_ok=True )
#         ofile = open( outfile, 'wb') 
#         pickle.dump(self , ofile , pickle.HIGHEST_PROTOCOL)
#         ofile.close()
#     
#         
# ###################################################################
#     def read_ts(self,ts_dir='.', verbose=True, name_filter='', add_key='', sites=[],lexclude=[], type = None, xyz=True):
# ###################################################################
#         """
#         Reads time series, trying to guess the format. Current time series format supported are: pos, kenv, mb_file, cats, txyz (pyacs), track (NEU format for high rate)
#         
#         :param ts_dir: directory of time series files
#         :param name_filter: string used to filter time series name '*name_filter*'
#         :param add_key: adds a string before site code
#         :param sites: list of site codes to be read. Any other will be discarded.
#         :param lexclude: list of sites to be excluded from reading
#         :param type: specifies the format of the time series files. Choose among ['pos', 'kenv', 'mb_file', 'cats', 'txyz', 'track' , 'pride','pck']
#         :param xyz: for pos files, reads the XYZ coordinates rather than dNEU. This is the default.
#         
#         :return: an Sgts instance.
#         """
#         # import
#         import numpy as np
# 
#         if verbose:
#             print('-- Reading directory: ',ts_dir)
# 
#         #######################################################################
#         def __pride_pos(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
#             
#             import glob
#             import os
#             
#             lGts={}
#             lfile = glob.glob(ts_dir+'/'+'pos_*'+name_filter)
#             lcode = []
#             for ifile in lfile:
#                 if not os.path.isfile( ifile ): continue
#                 code = ifile[-4:].upper()
#                 if code not in lcode:
#                     lcode.append( code )
#                 
#             for code in lcode:
#                 if verbose:
#                     print("-- Reading pride pos files for: %s " % code )
#                 lGts[code+add_key]=Gts(code=code)
#                 lGts[code+add_key].read_pride_pos(tsdir=ts_dir , verbose=verbose )
# 
#             return(lGts)
# 
#         #######################################################################
#         def __pride(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
#             import glob
#             import os
#             
#             lGts={}
#             lfile = glob.glob(ts_dir+'/'+'kin_*'+name_filter)
#             lcode = []
#             for ifile in lfile:
#                 if not os.path.isfile( ifile ): continue
#                 code = ifile[-4:].upper()
#                 if code not in lcode:
#                     lcode.append( code )
#                 
#             for code in lcode:
#                 if verbose:
#                     print("-- Reading pride files for: %s " % code )
#                 lGts[code+add_key]=Gts(code=code)
#                 lGts[code+add_key].read_pride(tsdir=ts_dir , verbose=verbose )
# 
#             return(lGts)
# 
#         #######################################################################
#         def __tdp(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
#             import glob
#             import os
#             
#             lGts={}
#             ldat=sorted(glob.glob(ts_dir+'/'+'*'+name_filter+'*.tdp'))
#             for fdat in ldat:
#                 if not os.path.isfile( fdat ): continue
#                 mb_file=fdat.split('/')[-1][0:4]
#                 code=mb_file[0:4]
#                 if code not in sites and sites !=[]:continue
#                 if verbose:
#                     print("-- Reading ", mb_file)
#                 lGts[code+add_key]=Gts(code=code)
#                 lGts[code+add_key].read_tdp(ifile=fdat,gmt=False)
# 
#             return(lGts)
# 
#         
#         #######################################################################
#         def __mb_files(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
#             import glob
#             import os
#             lGts={}
#             ldat=sorted(glob.glob(ts_dir+'/'+'*'+name_filter+'*.dat1'))
#             for fdat in ldat:
#                 if not os.path.isfile( fdat ): continue
#                 mb_file=fdat[0:-1]
#                 code=mb_file[-12:-8]
#                 if code not in sites and sites !=[]:continue
#                 if verbose:
#                     print("-- Reading ", mb_file)
#                 lGts[code+add_key]=Gts(code=code)
#                 lGts[code+add_key].read_mb_file(ifile=mb_file,gmt=False)
# 
#             return(lGts)
#         
# 
#         #######################################################################
#         def __kenv(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
#             import glob
#             import os
#             
#             lGts={}
#             lkenv=sorted(glob.glob(ts_dir+'/*.kenv'))
#             for fkenv in lkenv:
#                 if not os.path.isfile( fkenv ): continue
#                 code=fkenv[-20:-9]
#                 if verbose:print("-- Reading ", fkenv)
#                 lGts[code+add_key]=Gts(code=code)
#                 try:
#                     lGts[code+add_key].read_kenv(fkenv,date_type='jd')
#                 except:
#                     print("!!! Error. Could not read ",fkenv)
#                     del lGts[code+add_key]
#                     continue
#             return(lGts)
# 
#         #######################################################################
#         def __cats(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose):
#             import glob
#             import os
#             lGts={}
#             lcats=sorted(glob.glob(ts_dir+'/cats*.dat'))
#             for cats in lcats:
#                 if not os.path.isfile( cats ): continue
#                 code=cats.split('cats_')[-1][:4].upper()
#                 if verbose:print("-- Reading ", cats)
#                 lGts[code+add_key]=Gts(code=code)
#                 try:
#                     lGts[code+add_key].read_cats_file(ifile=cats,gmt=False)
#                 except:
#                     print("!!! Error. Could not read ",cats)
#                     del lGts[code+add_key]
#                     continue
#             return(lGts)
# 
#         #######################################################################
#         def __pos(ts_dir=ts_dir, name_filter=name_filter, add_key=add_key, verbose=verbose,xyz=True):
#             import glob
#             lGts={}
#             lpos=sorted(glob.glob(ts_dir+'/????*.pos'))
#             import os
#             for pos_w_path in lpos:
#                 if not os.path.isfile( pos_w_path): continue
#                 _drive, path_and_file = os.path.splitdrive(pos_w_path)
#                 _path, pos = os.path.split(path_and_file)
#                 code=pos[:4]
#                 if code not in sites and sites !=[]:continue
#                 if verbose:print("-- Reading ", pos)
#                 lGts[code+add_key]=Gts(code=code)
#                 lGts[code+add_key].read_pos(tsfile=pos_w_path,xyz=xyz)
#             return(lGts)
# 
#         #######################################################################
#         def __ts_xyz(ts_dir='.',name_filter=None, verbose=False):
#             
#             import glob
#             import os
#             
#             lGts={}
#             l_txyz=glob.glob(ts_dir+'/*txyz.dat')
#             
#                 
#             ##############################################
#             def __read_txyz(file_name,lGts,verbose=False):
# 
#                 """
#                 Reads a single txyz.dat
#                 """
#                 
#                 NP_STR_TXYZ=np.genfromtxt(file_name,dtype=str)
#                 NP_CODE=NP_STR_TXYZ[:,2]
#                 NP_DATE=np.array(NP_STR_TXYZ[:,1],dtype=float)
#                 NP_XYZ=np.array(NP_STR_TXYZ[:,4:],dtype=float)
#             
#                 for i in np.arange(NP_CODE.shape[0]):
#                     code = NP_CODE[i]
#                     
#                     if code in list(lGts.keys()):
#                         lGts[code].add_obs_xyz(NP_DATE[i],NP_XYZ[i],in_place=True,check=False)
#                     else:
#                         new_gts=Gts(code=code)
#                         new_gts.add_obs_xyz(NP_DATE[i],NP_XYZ[i],in_place=True,check=False)
#                         lGts[code]=new_gts
#             
#                 return(lGts)
#             
#             i=1
#             for txyz in l_txyz:
#                 if not os.path.isfile( txyz ): continue
#                 if verbose:
#                     print('-- Reading ',txyz,i,'/',len(l_txyz))
#                     i = i + 1
#                     
#                 lGts=__read_txyz(txyz, lGts,verbose=verbose)
#             
#             for code in sorted( lGts.keys() ):
#                 
#                 gts = lGts[code]
#                 
#                 if verbose:
#                     print('-- Generating NEU from XYZ for ',gts.code)
#                 gts.xyz2neu(corr=False)
#                 
#                 if not gts.cdata():
#                     sys.exit()
#                 
#             return(lGts)
#             
# 
#         #######################################################################
#         def __track_NEU(ts_dir='.', name_filter=name_filter, add_key=add_key, verbose=verbose):
#             import glob
#             lGts={}
#             ltrack=sorted(glob.glob(ts_dir+'/*.NEU.*'))
#             import os
#             for track_w_path in ltrack:
#                 if not os.path.isfile( track_w_path ): continue
#                 _drive, path_and_file = os.path.splitdrive(track_w_path)
#                 _path, track = os.path.split(path_and_file)
#                 code= track.split('.')[-2].upper()
#                 if code not in sites and sites !=[]:continue
#                 if verbose:print("-- Reading ", track)
#                 lGts[code+add_key]=Gts(code=code)
#                 try:
#                     lGts[code+add_key].read_track_NEU(tsfile=track_w_path)
#                 except:
#                     print("!!! Error. Could not read ",track)
#                     del lGts[code+add_key]
#                     continue
#             return(lGts)
# 
# 
#         # START
# 
#         lGts={}
#         #######################################################################
#         # PCK FORMAT
#         #######################################################################
#         # Added JMN 22/08/2019
#         # Sgts saved as a pickle dump
#         if type is None or type=='pck':
#             import glob
#             lpck = glob.glob( self.dir+'/*.pck' )
#             if lpck != []:
#                 pck = lpck[-1]
#                 # read Sgts from pickle
#                 import pickle
#                 with open( pck, "rb") as f:
#                     ts = pickle.load(f)
#                 f.close()
#                 
#                 for code in ts.lcode():
#                     lGts[code] = ts.__dict__[code]
#                 
#                 print("-- Sgts read from PYACS pck file: %s" % (pck) )
#             else:
#             # no pck file    
#                 if verbose:print("-- No PYACS pck file found")
# 
# 
#         #######################################################################
#         # PRIDE POS
#         #######################################################################
#         if type is None or type=='pride_pos':
#             lGts_pride_pos = __pride_pos(ts_dir=self.dir,verbose=verbose)
#             if lGts_pride_pos=={}:
#                 if verbose:print("-- No pride pos files found")
#             else:
#                 lGts=lGts_pride_pos
# 
#         #######################################################################
#         # PRIDE
#         #######################################################################
#         if type is None or type=='pride':
#             lGts_pride = __pride(ts_dir=self.dir,verbose=verbose)
#             if lGts_pride=={}:
#                 if verbose:print("-- No pride_files found")
#             else:
#                 lGts=lGts_pride
# 
#         #######################################################################
#         # MB_FILE
#         #######################################################################
#         if type is None or type=='mb_file':
#             lGts_mb_files=__mb_files(ts_dir=self.dir,verbose=verbose)
#             if lGts_mb_files=={}:
#                 if verbose:print("-- No mb_files found")
#             else:
#                 lGts=lGts_mb_files
#                 
#                 try:
#                     self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
#                 except:
#                     print("! Warning: Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files")
#         
#         #######################################################################
#         # TDP
#         #######################################################################
#         if type is None or type=='tdp':
#             lGts_mb_files=__tdp(ts_dir=self.dir,verbose=verbose)
#             if lGts_mb_files=={}:
#                 if verbose:print("-- No tdp_files found")
#             else:
#                 lGts=lGts_mb_files
#                 try:
#                     self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
#                 except:
#                     print("! Warning: Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files")
#         
#         #######################################################################
#         # KENV
#         #######################################################################
#         if type is None or type=='kenv':
#             lGts_kenv=__kenv(ts_dir=self.dir,verbose=verbose)
#             if lGts_kenv=={}:
#                 if verbose:print("-- No kenv file found")
#             else:lGts=lGts_kenv
#         
#         #######################################################################
#         # CATS
#         #######################################################################
#         if type is None or type=='cats':
#             lGts_cats=__cats(ts_dir=self.dir,verbose=verbose)
#             if lGts_cats=={}:
#                 if verbose:print("-- No cats file found")
#             else:
#                 lGts=lGts_cats
#                 try:
#                     self.read_gmt(verbose=verbose, gmt=self.dir+'/../stat/pyacs_void.gmt')
#                 except:
#                     print("! Warning: Could not read a gmt file to populate Gts attributes lon,lat,h & X0,Y0,Z0 for mb_files")
#         
#         
#         #######################################################################
#         # POS FORMAT
#         #######################################################################
#         if type is None or type=='pos':
#             lGts_pos=__pos(ts_dir=self.dir,verbose=verbose,xyz=xyz)
#             if lGts_pos=={}:
#                 if verbose:print("-- No Gamit/Globk pos file found")
#             else:lGts=lGts_pos
# 
#         #######################################################################
#         # TXYZ FORMAT
#         #######################################################################
#         if type is None or type=='txyz':
#             lGts_txyz=__ts_xyz(ts_dir=self.dir,verbose=verbose)
#             if lGts_txyz=={}:
#                 if verbose:print("-- No pyacs t_xyz file found")
#             else:lGts=lGts_txyz
# 
#         #######################################################################
#         # TRACK FORMAT
#         #######################################################################
#         if type is None or type=='track':
#             lGts_track=__track_NEU(ts_dir=self.dir,verbose=verbose)
#             if lGts_track=={}:
#                 if verbose:print("-- No Gamit/Globk track NEU file found")
#             else:lGts=lGts_track
# 
# 
#         for gts in list(lGts.values()):
#             if gts.code not in lexclude:
#                 print('-- adding ', gts.code)
#                 self.append(gts)
# 
#         if verbose:
#             print('-- read ',len(self.lcode() ),' time series in directory ',ts_dir)
#     
# 
# 
# ###################################################################
#     def read_gmt(self,gmt=True,verbose=False,vel=False):
# ###################################################################
#         """
#         Reads a gmt psvelo file and populates lon and lat attributes for each Gts of Sgts
#         
#         :param gmt: if True tries to read '/../stat/pyacs_void.gmt', if a string then it is the gmt file to be read. 
#         :param verbose: verbose mode
#         :param vel: boolean. If True, fills the .velocity attribute of every time series with the values read in the gmt file.
#         
#         :note: this method is always in place.
#         """
# 
#         if gmt==True:
#             gmt_file=self.dir+'/../stat/pyacs_void.gmt'
#         if  isinstance(gmt,str):
#             gmt_file=gmt
#         
#         if gmt != False:
#             if verbose:print("-- Importing modules Velocity_Field, Gpoint ",gmt_file)
#     
#             from pyacs.lib.vel_field import Velocity_Field
#             
#             if verbose:print("-- Reading ",gmt_file)
#             velf=Velocity_Field()
#             velf.read(file_name=gmt_file)
#             lsite=vel.lcode()
#             
#             for gts in self.lGts():
#                 if verbose:print("-- Processing ",gts.code)
#                 if gts.code in lsite:
#                     M=velf.site(gts.code)
#                     gts.lon=M.lon
#                     gts.lat=M.lat
#                     if vel:
#                         gts.velocity=[M.Ve,M.Vn,0.0,M.SVe,M.SVn,0.0]
# 
# 
# ###################################################################
#     def frame(self,frame=None,euler=None,w=None,verbose=False):
# ###################################################################
#         """
#         Rotates the time series according to an Euler pole.
#         User must provide either frame, euler or w.
#         
#         
#         :param frame: str, implemented values are 'soam','nas','nazca','inca','nas_wrt_soam','inca_wrt_soam'.
#         :param euler: Euler values provided either as a \
#         string 'euler_lon/euler_lat/euler_w', a list [euler_lon,euler_lat,euler_w] or \
#         a 1D numpy array np.array([euler_lon,euler_lat,euler_w])
#         :param w: rotation rate vector in rad/yr, provided either as a \
#         string 'wx/wy/wz', a list [wx,wy,wz] or \
#         a 1D numpy array np.array([wx,wy,wz])  
# 
#         :return: the new Sgts instance in new frame
# 
#         :ref: All values for frames are from Nocquet et al., Nat Geosc., 2014.
# 
#         """
# 
#         # import
#         import numpy as np
#         import pyacs.lib.euler
#         
#         # check arguments are OK
#         
#         if [frame,euler,w].count(None) != 2:
#             print('!!! ERROR: define either argument frame, euler or w ')
#             return(None)
#         
#         # Euler poles taken from pygvel_pole_info.py
#         lEuler={}
#         lEuler['soam']=[-132.21,-18.83,0.121]
#         lEuler['nas']=[-97.52,6.48,0.359]
#         lEuler['nazca']=[-94.4,61.0,0.57]
#         lEuler['inca']=[-103.729,-1.344,0.1659]
#         lEuler['nas_wrt_soam']=[-83.40,15.21,0.287]
#         lEuler['inca_wrt_soam']=[-63.76,22.47,0.092]
#         
#         # check frame case is OK
#         if ( frame not in list(lEuler.keys())) and ( frame is not None):
#             print("!!! ERROR: requested frame ",frame," not known")
#             print("!!! ERROR: available frames are: ", list(lEuler.keys()))
#             return(None)
#         
#         # initialize new gts
#         
#         New_Sgts=Sgts(read=False)
# 
#         # convert to Euler vector whatever the provided argument
#         
#         # case frame
#         if frame is not None:
#             euler_vector=np.array(lEuler[frame])
# 
#         # case w as rotation rate vector
#         if w != None:
# 
#             if ( isinstance(w,str) ) and '/' in w:
#                 w=np.array(list(map(float,w.split('/'))))
# 
#             if isinstance(w,list):
#                 w=np.array(w)
#             
#             if not isinstance(w,np.ndarray):
#                 print('!!! ERROR: argument w not understood: ',w)
#                 return(None) 
#             
#             euler_vector=np.array(pyacs.lib.euler.rot2euler([w[0],w[1],w[2]]))
# 
#         # case euler vector
#         if euler is not None:
# 
#             if ( isinstance(euler,str) ) and '/' in euler:
#                 euler=np.array(list(map(float,euler.split('/'))))
# 
#             if isinstance(euler,list):
#                 euler=np.array(euler)
#             
#             if not isinstance(euler,np.ndarray):
#                 print('!!! ERROR: argument euler not understood: ',euler)
#                 return(None) 
#             
#             euler_vector=np.array(euler)
#         
#         # converts the gts
#         
#         for gts in self.lGts():
#             if verbose:print("-- Processing ",gts.code)
#             try:
#                 new_gts=gts.remove_pole(euler_vector,pole_type='euler',in_place=False, verbose=verbose)
#             except (RuntimeError, TypeError, NameError):
#                 print("!!! Error processing ",gts.code)
#                 continue
#             if isinstance(new_gts,Gts):
#                 New_Sgts.append(new_gts)
#             else:
#                 print("!!! Error processing ",gts.code, "!!! No time series created.")
# 
#         return(New_Sgts)
# 
# ###################################################################
# ## copy Sgts
# ###################################################################
# 
#     def copy(self):
#         """
#         makes a (deep) copy of the Sgts object
#         
#         :return : a new Sgts instance
#         """
#         import copy
#         
#         new_Sgts=copy.deepcopy(self)
#         return(new_Sgts)
# 
# 
# ###################################################################
#     def read_gts_conf(self,gts_conf_file,verbose=False):
# ###################################################################
#         """
#         Reads a gts_conf_file
#         implemented commands in the file are:
#         #todo add_break      [site]  [date] # date is either [decyear] [doy year] [mday month year]
#         apply_offset   [site]  [offset_north,offset_east,offset_up] [date] # offset  applied is in mm, date is either [decyear] [doy year] [mday month year]
#         remove_offset   [site]  [offset_north,offset_east,offset_up] [date] # offset removed is in mm, date is either [decyear] [doy year] [mday month year]
#         #todo extract_periods [site] [date1,date2] # date is either [decyear] [doy year] [mday month year]
#         #todo exclude_periods [site] [date1,date2] # date is either [decyear] [doy year] [mday month year]
#         #todo remove_day    [site] [date]
#         """
#         # import
#         import numpy as np
#         from pyacs.lib.astrotime import guess_date
# 
#         New_Sgts=self.copy()
# 
#         conf=open(gts_conf_file,'r')
# 
#         H_apply_offsets={}
#         H_offsets_dates={}
#         
#         for line in conf:
#             if len(line)<5 or line[0]=='#':continue
#             lline=line.split()
#             # apply_offset
#             if lline[0] in ['apply_offset','remove_offset']:
#                 (code,sdn,sde,sdu,sdate)=lline[1:]
#                 if lline[0] == 'apply_offset':
#                     dn=float(sdn)/1000.0
#                     de=float(sde)/1000.0
#                     du=float(sdu)/1000.0
# 
#                 if lline[0] == 'remove_offset':
#                     dn = -float(sdn)/1000.0
#                     de = -float(sde)/1000.0
#                     du = -float(sdu)/1000.0
# 
#                 
#                 date=guess_date(sdate)
#                 if code in H_apply_offsets:
#                     if verbose:print('-- adding offset to ',code,date,dn,de,du)
#                     H_apply_offsets[code]=np.vstack((H_apply_offsets[code],np.array([date,dn,de,du])))
#                     H_offsets_dates[code].append(date)
#                 else:
#                     print('-- adding new offset to ',code,date,dn,de,du)
#                     H_apply_offsets[code]=np.array([[date,dn,de,du]])
#                     H_offsets_dates[code]=[date]
#         
#         # apply commands
#         for gts in self.lGts():
#             if verbose:print("-- Processing ",gts.code)
#             if gts.code in H_apply_offsets:
#                 #gts.offsets_dates+=H_offsets_dates[gts.code]
#                 try:
#                     New_Sgts.__dict__[gts.code]=gts.apply_offsets(H_apply_offsets[gts.code],verbose=verbose)
#                     if verbose:
#                         print("-- applying offset ", H_apply_offsets[gts.code])
#                 except (RuntimeError, TypeError, NameError):
#                     print("!!! Error applying offset to ",gts.code)
#                     continue
# #            else:
# #                new_gts=Gts.copy(gts)
# #            if isinstance(new_gts,Gts):
# #                New_Sgts.append(new_gts)
# #            else:
# #                print "Error processing ",gts.code, "!!! No time series created."
# 
#         return(New_Sgts)
#         
# 
#         
# ###################################################################
#     def stat_site(self,lsite=[],lsite_file=None,verbose=False, save_dir=None):
# ###################################################################
#         """
#         basic statistics about individual time series
# 
#         :param lsite: list of selected sites for statistics 
#         :param lsite_file: list of selected sites for statistics provided as a file
#         :param verbose: verbose mode
#         :param save_dir: directory where statistics files will be written 
#         
#         """
#         # import
#         import numpy as np
#         import os
#         
#         # Case lsite_file option
#         
#         if lsite_file:
#             f=open(lsite_file,'r')
#             for line in f:
#                 lline=line.split()
#                 if (len(lline[0])!=4):continue
#                 lsite.append(line.split()[0])
#             f.close()
#         
#         
#         # save_dir case
# 
#         if save_dir is not None:
# 
#             if os.path.isfile(save_dir):
#                 print('!!! ',save_dir,' file already exists')
#                 sys.exit()
#             
#             else:
#                 if not os.path.isdir(save_dir):
#                     try:
#                         os.mkdir(save_dir)
#                     except:
#                         print('!!! ERROR: Could not create ' , save_dir)
#                         sys.exit()
# 
#             info_file=open(save_dir+'/pyacs.info','w')
# 
#         
#         # header & print formats
#         
#         header   = "site     long.      lat.     v_e     v_n     v_u    sv_e    sv_n    sv_u   #obs #camp   s_year   e_year d_year     se     sn     su   wrms_e   wrms_n   wrms_u"
#         sep_line = "--------------------------------------------------------------------------------------------------------------------------------------------------------------"
#         fmt      = "%4s %9.4lf  %8.4lf %7.2lf %7.2lf %7.2lf %7.2lf %7.2lf %7.2lf %6d %5d %8.3lf %8.3lf  %5.2lf %6.2lf %6.2lf %6.2lf %8.2lf %8.2lf %8.2lf"
#         fmt_vna  = "%4s %9.4lf  %8.4lf     N/A     N/A     N/A     N/A     N/A     N/A %6d %5d %8.3lf %8.3lf  %5.2lf %6.2lf %6.2lf %6.2lf %8.2lf %8.2lf %8.2lf"
#         fmt_na   = "%4s %9.4lf  %8.4lf     N/A     N/A     N/A     N/A     N/A     N/A %6d %5d %8.3lf %8.3lf  %5.2lf    N/A    N/A    N/A      N/A      N/A      N/A"
#         
#         if save_dir:
#             info_file.write(header+"\n")
#             info_file.write(sep_line+"\n")
#         
#         print(header)
#         print(sep_line)
#         # loop on Gts
#         
#         for gts in self.lGts():
#             
#             # if no data move to next gts
#             if gts.data is None:
#                 continue
#             
#             if lsite != [] and gts.code not in lsite:continue
#             code=gts.code
#             lon=gts.lon
#             lat=gts.lat
#             n_obs=gts.data.shape[0]
#             start_year=np.min(gts.data[:,0])
#             end_year=np.max(gts.data[:,0])
#             delta_year=end_year-start_year
#             n_camp=np.unique(list(map(int,list(gts.data[:,0])))).shape[0]
#             
# #            if delta_year>2 and n_obs>2 and n_camp>1:
#             
#             # detrend
#             new_gts=gts.detrend_median( auto=True )
# 
#             if new_gts is not None:
#                 
#                 [wrms_n,wrms_e,wrms_u]=new_gts.data[:,1:4].std(axis=0)*1000.0
#                 [s_n   ,   s_e,s_u   ]= np.median( np.fabs( np.diff( new_gts.data[:,1:4]*1000.0 , axis=0 ) ) , axis=0 )
#                 [vn,ve,vu,svn,sve,svu]=new_gts.velocity[:6]*1000.0
#                 
#                 fmt_str = (fmt % (code,lon,lat, ve,vn,vu, sve, svn, svu, n_obs,n_camp, start_year,end_year,delta_year,s_n,s_e,s_u,wrms_e,wrms_n,wrms_u))
#                 print(fmt_str)
#                 
#                 if save_dir is not None:
#                     info_file.write("%s\n" % (fmt_str))
#                     
#             else:
#                 if gts.data.shape[0] >= 2:
#                     # we can still get wrms and daily scatter
#                     [wrms_n,wrms_e,wrms_u]= gts.data[:,1:4].std(axis=0)*1000.0
#                     [s_n   ,   s_e,s_u   ]= np.median( np.fabs( np.diff( gts.data[:,1:4]*1000.0 , axis=0 ) ) , axis=0 )
# 
#                     fmt_str = (fmt_vna % (code,lon,lat, n_obs,n_camp, start_year,end_year,delta_year,s_n,s_e,s_u,wrms_e,wrms_n,wrms_u))
#                     print(fmt_str)
# 
#                     if save_dir is not None:
#                         info_file.write("%s\n" % (fmt_str))
#                 
#                 else:
#                     fmt_str = (fmt_na % (code,lon,lat, n_obs,n_camp, start_year,end_year,delta_year))
#                     print(fmt_str)
# 
#                     if save_dir is not None:
#                         info_file.write("%s\n" % (fmt_str))
#                     
# 
#         # write results
#         
#         if save_dir is not None:
# 
#             info_file.close()
#             info = np.genfromtxt( save_dir+'/pyacs.info' , skip_header=2,dtype=str )
# 
# 
#             header_info_ds = "site    long.       lat.\n"\
#                            + "------------------------"
#             np.savetxt(save_dir+'/pyacs_lphn.dat',info[:,:3] , fmt = "%4s %10s %10s" , header=header_info_ds )
#             
#             info_ds = info[ np.where(info[:,14]!='N/A') ]
#             header_info_ds = "site daily_scatter_north (mm) daily_scatter_east (mm) daily_scatter_up (mm)\n"\
#                            + "---------------------------------------------------------------------------"
#             np.savetxt(save_dir+'/pyacs_daily_scatter.dat',info_ds[:, (0,14,15,16) ] , fmt = "%4s %20s %20s %20s" , header=header_info_ds)
#             
#             info_ds = info[ np.where(info[:,17]!='N/A') ]
#             header_info_ds = "site    wrms_north       wrms_east    wrms_up (mm) \n"\
#                            + "---------------------------------------------------"
#             np.savetxt(save_dir+'/pyacs_wrms.dat',info_ds[:, (0,17,18,19) ] , fmt = "%4s %15s %15s %15s" , header=header_info_ds)
#             
#             np.savetxt(save_dir+'/pyacs_void.dat',info[:, (1,2,0) ] , fmt = "%10s %10s   0.00   0.00   0.00   0.00    0.00   %s")
# 
#             info_ds = info[ np.where(info[:,3]!='N/A') ]
#             header_info_ds = "   long.       lat.         ve         vn        sve        svn    sven    site\n"\
#                            + "------------------------------------------------------------------------------------"
#             np.savetxt(save_dir+'/pyacs_vel.dat',info_ds[:, (1,2,3,4,6,7,0) ] , fmt = "%10s %10s %10s %10s %10s %10s    0.00    %s", header=header_info_ds)
# 
#             header_info_ds = "   long.       lat.        --       v_up       ---   sigma_v_up      ---    site\n"\
#                            + "---------------------------------------------------------------------------------"
#             np.savetxt(save_dir+'/pyacs_vel_up.dat',info_ds[:, (1,2,5,8,0) ] , fmt = "%10s %10s      0.00 %10s      0.00 %12s     0.00    %s" , header=header_info_ds)
#         
#         
#         lcode=[]
#         for gts in self.lGts():lcode.append(gts.code)
#         for site in lsite:
#             if site not in lcode:print("!!! ",site," in ",lsite_file," and not in ts")
#         
# ###################################################################
#     def add_offsets_dates(self,dates,verbose=False):
# ###################################################################
#         """
#         add_offsets_dates to every Gts in current Sgts
#         """
#         New_Sgts=Sgts(read=False)
#         for gts in self.lGts():
#             if verbose:print("-- Processing ",gts.code)
#             try:
#                 new_gts=gts
#                 new_gts.offsets_dates=dates
#             except (RuntimeError, TypeError, NameError):
#                 print("!!! Error processing ",gts.code)
#                 continue
#             if isinstance(new_gts,Gts):
#                 New_Sgts.append(new_gts)
#             else:
#                 print("!!! Error processing ",gts.code, "!!! No time series created.")
#         return(New_Sgts)
# 
# ###################################################################
#     def read_soln(self,soln,verbose=True):
# ###################################################################
#         """
#         read a IGS soln file and add an offsets_dates for any P change in soln.
#         
#         :param soln: soln.snx IGS file
#         
#         :note: the method is in place
#         """
#     
#         from pyacs.sol.discontinuity import Discontinuities
#         
#         fsoln=Discontinuities()
#         if verbose:print("-- Reading ",soln)
#         fsoln.read_igs_discontinuity(soln)
#     
#         for gts in self.lGts():
#             if verbose:print("-- Adding offsets to ",gts.code)
#             if gts.code in fsoln.lsite_offsets_dates:
#                 gts.offsets_dates=fsoln.lsite_offsets_dates[gts.code]
#     
#         return()
#         
# ###################################################################
# ##  DEFAULT PLOT FOR ALL TIMES SERIES
# ###################################################################
# 
#     def plot_all(self,odir='.',save=True,show=False,verbose=True,center=False,date=[],):
#         """
#         Creates a default plot for all time series provided as a dictionary
#         save=True creates png files (default True)
#         show=True show time series in matplotlib graphic window (default False)
#         
#         :param odir: output directory
#         :save : boolean. True (default) will create png files
#         :show : boolean. show time series in matplotlib graphic window (default False)
#         :verbose: boolean. Verbose mode
#         :center: boolean. If True, the time series will be centered around 0. Default is False.
#         :date: list [start_date,end_date] in decimal year. Default is [] and will adapt the date scale according to the data
#         
#         """
#     
#         for gts in self.lGts():
#             if verbose:print("-- Plotting ",gts.code)
#             if save:
#                 png=gts.code+'.png'
#             gts.plot(date=date,save=png,show=show,center=center)
#     
#         return()
# 
# ###################################################################
# ##  save velocity
# ###################################################################
#     def save_velocity(self,vel_file='../stat/vel'):
#         """
#         save horizontal and up velocities
#         """
#         h_vel=vel_file+'_en.gmt'
#         u_vel=vel_file+'_up.gmt'
#         for gts in self.lGts():
#             try:
#                 gts.save_velocity(h_vel)
#             except:
#                 print("!!! Could not save ",gts.code," in ",h_vel)
#             try:
#                 gts.save_velocity(u_vel,up=True)
#             except:
#                 print("!!! Could not save ",gts.code," in ",u_vel)
#         return()
# 
# ###################################################################
# ##  LAZY PYACS
# ###################################################################
# 
# #     def lazy_pyacs(self,wdir='../lazy',save=True,verbose=True):
# #         """
# #         Run lazy_pyacs on all time series
# #         """
# #         import os.path
# #         if not os.path.isdir(wdir):
# #             try:
# #                 os.mkdir(wdir)
# #                 if verbose:print("=> Writing lazy_pyacs results into dir ",wdir)
# #             except:
# #                 print("!!! Could not create directory: ",wdir)
# #                 
# #         lazy_Gts=Sgts(read=False)
# #         
# #         for gts in self.lGts():
# #             if verbose:print("=> Lazy pyacs on ",gts.code)
# #             try:
# #                 my_lazy_ts=gts.lazy_pyacs()
# #                 lazy_Gts.append(my_lazy_ts.remove_outliers().detrend())
# #                 my_lazy_ts.remove_outliers().detrend().write_pos(add_key='lazy',ts_dir=wdir)
# #             except:
# #                 print("!!! Could not complete lazy_pyacs analysis")
# #         
# #         lazy_Gts.save_velocity()
# #     
# #         return()
# 
# ###################################################################
#     def sel_rectangle(self,bounds,verbose=True):
# ###################################################################
#         """
#         selects the time series for sites within a rectangles
#         
#         :param bounds: [lon_min,lon_max,lat_min,lat_max]
#         :pram verbose: verbose mode
#         
#         :return: a new Sgts instance
#         """
#         [lon_min,lon_max,lat_min,lat_max]=bounds
# 
#         new_Sgts=Sgts(read=False)
#         for gts in self.lGts():
#             current_lon=gts.lon
#             current_lat=gts.lat
#             
#             if current_lon>=lon_min and current_lon<=lon_max<=lon_max and current_lat>= lat_min and current_lat<=lat_max:
#                 if verbose:print("-- ",gts.code," selected")
#                 new_Sgts.append(gts)
#         return(new_Sgts)
#         
# 
# ###################################################################
#     def to_displacement(self,verbose=True,base_name='vel',wdir='.',up=False):
# ###################################################################
#         """
#         print displacements every dates as gmt psvelo files 
#         """
#         
#         # import
#         
#         import pyacs.lib.astrotime as AT
#         import pyacs.lib.vel_field as Velocity_Field
#         import pyacs.lib.gmtpoint as GMT_Point
#         import numpy as np
#         
#         # get the list of dates
#         
#         ldate=[]
#         for gts in self.lGts():
#             ldate=ldate+gts.data[:,0].tolist()
# 
#         np_dates=np.array(sorted(list(set(ldate))))
#         
#         if verbose: print("=> Found ",np_dates.shape[0]," dates")
# 
#         # creates the gmt files
#         
#         date_ref=np_dates[0]
#         
#         for i in np.arange(np_dates.shape[0]):
#             date=np_dates[i]
# #        for date in sorted(ldate):
#             file_name=wdir+'/'+("%04d_" % i)+base_name+'.gmt'
#             vel=Velocity_Field.Velocity_Field(file_name=file_name)
#             for gts in self.lGts():
#                 M=GMT_Point.GMT_Point(code=gts.code)
#                 
#                 tmp_gts=gts.extract_dates([date])
#                 
#                 if tmp_gts.data is not None:
#                     if tmp_gts.data.shape[0]==1:
#                         M.lon=gts.lon
#                         M.lat=gts.lat
#                         if up:
#                             M.Ve=0.0
#                             M.Vn=tmp_gts.data[0,3]*1.E3
#                             M.SVe=0.0
#                             M.SVn=tmp_gts.data[0,6]*1.E3
#                             M.SVen=0.0
#                         else:
#                             M.Ve=tmp_gts.data[0,2]*1.E3
#                             M.Vn=tmp_gts.data[0,1]*1.E3
#                             M.SVe=tmp_gts.data[0,5]*1.E3
#                             M.SVn=tmp_gts.data[0,4]*1.E3
#                             M.SVen=0.0
#                         
#                         vel.add_point(M)
#             
#             (mday,month,ut)=AT.decyear2cal(date)
#             (noday,ut)=AT.decyear2dayno(date)    
#             date_info=("step #%04d date %10.5lf %02d-%02d-%04d-%.1lf doy %03d | days since reference date +%4.1lf " % \
#                        (i,date,mday,month,int(date),ut,noday,(date-date_ref)*365.25))
# 
#             vel.write(comment=date_info)
# 
# ###################################################################
# ##  MEDVEL
# ###################################################################
# 
#     def medvel(self,outdir=None,verbose=False):
#         """
#         
#         Automatic velocity estimates using median estimator.
#         The code is adapted from the MIDAS approach (Blewitt et al., 2016).
#         
#         medvel fills the velocity attribute of every Gts from the current Sgts instance.
#         
#         returns the modified Sgts instance
#         Optionally, if outdir option is provided, writes the results in outdir
#         
#         :param: outdir: output directory, default None
#         :param: verbose: boolean, verbose mode
#         :param: warning: output warning file
# 
#         :reference: Blewitt, G., Kreemer, C., Hammond, W. C., & Gazeaux, J. (2016). MIDAS robust trend estimator for accurate GPS station velocities without step detection. Journal of Geophysical Research: Solid Earth, 121(3), 2054-2068.
#         """
# 
#         import os
#         
#         # 
#         
#         new_sgts = Sgts(read=False)
#         
#         # create output directory
#         
#         if outdir is not None:
#             
#             if os.path.isdir(outdir):
#                 print('!!! ',outdir,' directory already exists')
#                 sys.exit()
#             
#             elif os.path.isfile(outdir):
#                 print('!!! ',outdir,' file already exists')
#                 sys.exit()
#             
#             else:
#                 if verbose:
#                     print('-- creating ' , outdir )
#                 
#                 os.mkdir(outdir)
#             
#             # output file names
#             
#             out_cgps = outdir+'/vel_cgps.dat'
#             out_sgps = outdir+'/vel_sgps.dat'
#     
#             out_cgps_up = outdir+'/vel_cgps_up.dat'
#             out_sgps_up = outdir+'/vel_sgps_up.dat'
#             
#             # open warning file
#             
#             fwarning = open(outdir+'/warning.dat','w+')
#         
#         # start loop on sites
#         
#         lcode = self.lcode()
#         
#         for site in lcode:
# 
#             
#             if verbose:
#                 print('-- Processing ', site )
# 
#             # check whether there is at least one year of data
#             
#             if ( self.__dict__[site].data[-1,0] - self.__dict__[site].data[0,0] ) < 1.0 :
#                 if verbose:
#                     print("-- Less than one year of data for site: %s" % site)
#                 if outdir is not None:
#                     fwarning.write("-- Less than one year of data for site: %s\n" % site)
#                 continue
# 
#             # check whether there are at least three data
#             
#             if ( self.__dict__[site].data.shape[0] ) < 3 :
#                 if verbose:
#                     print("-- Less than 3 data for site: %s" % site)
#                 if outdir is not None:
#                     fwarning.write("-- Less than 3 data for site: %s\n" % site)
# #                continue
# 
#             # 
# 
#             detrended = self.__dict__[site].detrend_median( auto=True )
# 
#             if ( detrended.velocity is None ):
#                 print( "-- Problem in detrend_median for site : %s " % site )
#                 continue
# 
#             if  ( outdir is not None ):
# 
#                 detrended.save_velocity(out_cgps,verbose=verbose)
#                 detrended.save_velocity(out_cgps_up,verbose=verbose,up=True)
#                 new_sgts.append( detrended )
#         
#         return( new_sgts )       
# 
#                     
# ###################################################################
# ##  SAME_SITE
# ###################################################################
# 
#     def same_site(self,dc=10, in_place=True, verbose=False):
#         """
#         
#         Check that all gts in the current Sgts are actually the same site. If a given time series is
#         found to be of two separate sites, then a new gts is added to the return Sgts instance.
# 
#         param dc: critical distance to decide to split the time series
#         param in_place: if True modify current Sgts, False retuen a new Sgts
#         param verbose: verbose mode
#         
#         return: a new Sgts instance
#         """
#         
#         # import
#         import numpy as np
# 
#         if not in_place:
#             new_Sgts = Sgts(read=False)
# 
#         # start loop on sites
#         
#         lcode = self.lcode()
#         
#         for site in lcode:
#             
#             if verbose:
#                 print('-- Processing ', site )
# 
#             my_ts = self.__dict__[site].copy()
# 
#             if my_ts.data_xyz is not None:
#                 data=my_ts.data_xyz[:,1:4]
#                 ddata=np.copy(my_ts.data_xyz[:,1:4])
#             
#             else:
#                 # if no data_xyz go to next gts
#                 print("!!! WARNING: data_xyz attribute required for method same_site and not found gts %s" % (site))
#             
#             # ensure median calculation
#             if np.mod(data.shape[0],2)==0:
#                 # duplicates the last date
#                 ddata = np.vstack((ddata,ddata[-1,:]))
#             
#             median = np.median(ddata,axis=0)
#             dist_data= np.sqrt( np.sum( (data-median)**2,axis=1) )
#             
#             lindex = np.where(dist_data > dc*1.E3 )
#             
#             # case gts needs to be split
#             if len( lindex[0] ) > 0 :
#                 # create a new code
#                 
#                 new_code = my_ts.code[:3]+'_'
#                 if new_code in self.lcode():
#                     import sys
#                     print("!!! ERROR: try to create a new gts with code %s and it already exists." % (new_code))
#                     new_code = my_ts.code[:2]+'__'
#                     
#                 print("-- time series for site %s appears to include different sites because there are coordinates at %d dates %.1lf km from the median position" % ( site, len(lindex) , np.max( ddata )*1.E-3 ) ) 
#                 print("-- %s time series will be split into code %s and code %s" % (site,site,new_code) )
#                 
#                 # create a new gts
#                 new_gts = Gts(code=new_code,data_xyz=np.copy(my_ts.data_xyz[lindex]))
#                 new_gts.xyz2neu(corr=True)
# 
#                 # remove the line from my_ts                
# 
#                 my_ts.data_xyz = np.delete( my_ts.data_xyz , lindex , axis=0 )
# 
#                 my_ts.xyz2neu(corr=True)
#                 
#                 # update the ouput
#                 
#                 if in_place:
#                     self.append(new_gts)
#                 else:
#                     new_Sgts.append( new_gts )
#             
# 
#             if in_place:
#                 self.__dict__[site] = my_ts
#             else:
#                 new_Sgts.append( my_ts )
# 
#         if in_place:
#             return self
#         else:
#             return new_Sgts
#             
#         
#         
#         
#         
#         
#         
#         
        