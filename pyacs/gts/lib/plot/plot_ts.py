###############################################################################
def plot_ts( ax, np_date , data, yerr, \
            loutliers=[], \
            yaxis=None, \
            error_scale=1.0,\
            **H_kwargs):
###############################################################################

    # import
    import pyacs.lib.astrotime
    import numpy as np

    if yaxis:
        (ymin,ymax)=yaxis
        ax.set_ylim( ymin , ymax )
    

    # plot data
    ax.errorbar(np_date,data, yerr=yerr*error_scale, **H_kwargs)

    # plot outliers
    # not implemented yet
#     if loutliers:
#         ts_outliers = __window_ts(data[loutliers,:],xmin=xmin,xmax=xmax)
#         if ts_outliers==[]:return()
#         
#         ts_outliers[:,1]=ts_outliers[:,1]-yaxis_shift
#         
#         plt.errorbar(ts_outliers[:,0],ts_outliers[:,1], yerr=ts_outliers[:,2]*error_scale, fmt='o', markersize=H_kwargs['markersize'], color='r',ecolor='grey',capsize=0,linewidth=1.0)
    
    ax.grid(True)
    

    return
