from setuptools import setup

###############################################################################
# VERSION HISTORY
# 0.64.0 on 20200910 : write pck now use DEFAULT_PROTOCOL instead of HIGHEST for backward compatibility
# 0.63.9 on 20200728 : working version, plwf & prox-tv removed from requirement in setup.py
# 0.63.9 on 20200508 : working version
# 0.63.8 on 20200325 : refactoring version glinalg
# 0.63.7 on 20200212 : add common mode testing status
# 0.63.6 on 20200207 : add make_normal_system in pyacs.lib.glinalg
# 0.63.5 on 20200104 : Sgts reorganized
# 0.63.4 on 20200104 : handling of .data and .data_xyz is hopefully more consistent
# 0.63.3 on 20191217 : some refactoring. gts.read constructor. Should be compatible with jupyter notebook plotting_time_series.ipynb
#                      - turns display off by default
# 0.63.2 on 20191218 : gts.plot improvement.
#                      - turns display off by default
# 0.63.1 on 20191010 : add seconds conversion in pyacs.lib.astrotime
# 0.63.0 on 20190819 : add option ref_only_all in pyacs.sol.read_conf, pyacs.sol_sinex, pyacs_make_time_series.py
#                      fix non optimized tsxyz approach
#                      pck format for Sgts added ; also added in paycs_make_time_series.py
#                      pyacs.gts.plot: add info option , superimposed can have several gts 
# 0.62.9 on 20190814 : refactoring of gts.lib.primitive. Fix bug handling .data_xyz
# 0.62.8 on 20190805 : added pyacs.lib.astrotime & coordinates import in pyacs.py. Added pyacs.lib.utils. 
#                      Bug corrected in pyacs.lib.shapefile when wrting polygons (due to syntax change in pyshp).
# 0.62.7 on 20190607 : added gts.primitive.substract_ts_daily
# 0.62.6 on 20190604 : added spherical_baseline_length_rate in pyacs.lib.coordinates
# 0.62.5 on 20190520 : added scripts/pyacs_gvel_strain.py
# 0.62.4 on 20190514 : Gts.lib.format.py re-organized in Gts.lib.format
# 0.62.3 on 20190503 : reorganize filter for gts. Added total variation filters and piecewise linear fit.
# 0.62.2 on 20190411 : in gts.extract_periods if data_xyz present then data is rebuilt
#                    : and X0,Y0,Z0,lon,lat,h attributes are updated. Bug in data_xyz extraction corrected.
# 0.62.1 on 20190409 : added pyacs_gvel_pole.py script
# 0.62.0 on 20190409 : small potential bugs in astrotime corrected. Enforcing output type to int for mday, month, noday
# 0.62.0 on 20190401 : pygeca separated from pyacs
# 0.61.6 on 20190312 : script for combining psvelo velocity field added. Experimental.
# 0.61.5 on 20190124 : minor bug in gts.remove_outliers when max index of outliers exceeds the length of the new time series 
# 0.61.5 on 20190124 : minor bug in paycs_make_time_series.py where the correlation coefficients were not properly filled by rotating the XYZ cov matrix 
# 0.61.5 on 20190123 : small change in Sgts.plot to handle the save option
# 0.61.5 on 20190122 : change in Sgts.medvel. Now return the Sgts instance with velocity populated
# 0.61.5 on 20181227 : small bug in gts.plot which was preventing to see the first point with date_unit='cal' option ; also dates are plot at their true time rather at the day at 00:00
# 0.61.5 on 20181220 : working version for pyeq_kinematic_inversion.py python 3.7
# 0.61.4 on 20181217 : first Exception handling for gts. Still in test
# 0.61.3 on 20181211 : scaling_laws and magnitude re-introduced, pole estimates, work on median & medvel in Sgts
# 0.61.3 on 20181120 : working version for Green's function cleaning
# 0.61.2 on 20181112 : new calc_pole method in vel_field
# 0.61.2 on 20181112 : bug in Sgts.stat_site where ve and vn had been swaped
# 0.61.2 on 20181112 : bug in gts.exclude_periods corrected
# 0.61.2 on 20181005 : non linear trajectory model added to gts
# 0.61.1.20181004 : Pure python Okada's formulas implementation ; bug corrected in option periods in gts.detrend_median() 
# 0.61.1.20180927 : Major change in pyacs_make_time_series.py that now allows coordinates uncertainties to be calculated through projection
# 0.60.1.20180902: first python 3.6 working bundle with pyacs_make_time_series and the core library ok. 
# 0.60: - for version > 0.60, PYACS is being progressively migrated to python 3.6 & 3.7, only few features 0.51/python2.7 available for now 
# 0.51: - working code with pygeca
# 0.50: - added observation reweighting options. InSAR capability tested. 
# 0.49: - added insar capability in pyeq_static_inversion.py (quick made from Amazonia, after World Cup final, requires more carefully result printing) 
# 0.48: - add pyeq_model_to_disp.py
# 0.47: - adding obs-model displacement psvelo gmt files in pyeq_kinematics.py
#     : - occasional bug in apply_offset corrected (> dev3)
#     : - bug in coordinates.xyz_spherical_distance (arcos missing) 
# 0.46: - adding pyeq_scaling_laws.py
#     : - adding pyacs/lib/GMT.py for backwards compatibility (> dev2)
# 0.45: - corrected sign error for up Green function using tde in 
#     : - detrend_median added to Gts.model
#     : - bug gts.file -> gts.ifile (for > dev2) 
#     : - automatic shapefiles for cumulative slip and plots of time series in pyeq_kinematic.py/print_log (dev >3)
# 0.44: - datetime_from_calarray added to pyacs.lib.astrotime
#     : - bug corrected in uts2hmsmicros for seconds calculation
#     : - read_tdp added to gts.lib.format
#     : - added itermax=10 option in Gts.lazy
#     : - added force_day to gts to force daily solution at 12:00
#     : - change to gts.format.write_pos to force round at the second
# 0.43: - prototype of new detect_steps methods
# 0.42: - new lazy command with offset detection using a median filter implemented
#        - bug correction in plot when use of date and date_unit='days'
# 0.41: - lazy_pyacs seems to be operational again. Patch in pyeq_parametrize_curve_surface_triangles.py and coordinates for GMT5 compatibility
# 0.40: - geometry and green generation scripts included, as well as the static inversion -not fully operational yet.
# 0.39: - new gts method in Sgts
# 0.38: - this version has been distributed to Vergnolle/DeChabalier/Tissandier/
# 0.35: - slip time dependent modeling operational again
#       - pyacs_qgis_model2polygon.py script added to distribution
###############################################################################
# TODO
# 04/10/2018 : detrend_median_seasonal: periods option to be implemented 
#        think how seasonal terms are then removed for the part not used in parametric seasonal signal 
###############################################################################


setup(name='pyacs',
      version='0.64.0',
      description='PYACS: Geodetic analysis and modeling tools for Tectonics',
      long_description='Geodetic analysis and modeling tools for Tectonics.',
 
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Geodesy :: Geophysics',
      ],
      keywords='geodesy GPS earthquake elastic dislocation tectonics time series',
      url='',
      author='Jean-Mathieu Nocquet (Geoazur, IRD, CNRS, OCA, Cote d Azur University, France)',
      author_email='nocquet@geoazur.unice.fr',
      license='NA',
      packages=['pyacs',\
                'pyacs.lib',\
                # gts
                'pyacs.gts','pyacs.gts.lib', \
                'pyacs.gts.lib.non_linear_model' , \
                'pyacs.gts.lib.filters' , \
                'pyacs.gts.lib.format',\
                'pyacs.gts.lib.primitive',\
                'pyacs.gts.lib.plot',\
                # sol & sinex
                'pyacs.sol','pyacs.sinex',\
                # Sgts
                'pyacs.gts.Sgts_methods',\
                # refactoring 0.63.8
                'pyacs.glinalg'
                ],

       scripts=[\
                # IPYACS
                'pyacs/scripts/ipyacs.py',
                # PYACS MAKE TIME SERIES
                'pyacs/scripts/pyacs_make_time_series.py', \
#                # QGIS TRANSLATERS
                'pyacs/scripts/pyacs_qgis_psvelo_2_shapefile.py',\
#                # PYACS GAMIT
#                # PYACS GVEL POLE & STRAIN
                'pyacs/scripts/pyacs_gvel_pole.py',\
                'pyacs/scripts/pyacs_gvel_estimate_pole.py',\
                'pyacs/scripts/pyacs_gvel_prep_comb.py',\
                'pyacs/scripts/pyacs_gvel_comb.py',\
                'pyacs/scripts/pyacs_gvel_strain.py',\
#                
                ],
      install_requires=['ipython',   \
                        'numpy',     \
                        'scipy',     \
                        'matplotlib',\
                        'argparse',  \
                        'pyshp',     \
                        'pyaml',     \
                        'ansicolors', \
                        'pyshp>=2.0.1',
#                        'pwlf' , \
#                        'prox_tv'
                        ],
#                        'hdf5',  \
#                        'netCDF4'],
      zip_safe=False,

      test_suite='nose.collector',
      tests_require=['nose'],
)
