def show_map( self , bounds = None , highlight=[]):
    """
    Show a simple map of site location
    """
    import geopandas
    import matplotlib.pyplot as plt
    import numpy as np

    #fig = plt.figure()
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # get coordinates
    llong = []
    llat  = []
    lcode = []
    for code in sorted( self.lcode() ):
        lcode.append( code )
        llong.append( self.__dict__[code].lon )
        llat.append( self.__dict__[code].lat )

    # get bounds if not provided
    if bounds is None:
        lon_min = np.min( llong )
        lon_max = np.max( llong )
        lat_min = np.min( llat )
        lat_max = np.max( llat )
    else:
        [lon_min,lon_max,lat_min,lat_max] = bounds 

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world.plot()
    plt.xlim(lon_min,lon_max)
    plt.ylim(lat_min,lat_max)
    plt.plot(llong,llat,'ro',markersize=3)
  
    for i, txt in enumerate(lcode):
        plt.annotate(txt, (llong[i], llat[i]), fontsize=10)

    for code in highlight:
        plt.annotate(code, (self.__dict__[code].lon, self.__dict__[code].lat), fontsize=10 , color='r')
        plt.plot( [self.__dict__[code].lon], [self.__dict__[code].lat], 'ws',markersize=5)

    
    # show
    plt.show()
