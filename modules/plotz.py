""" Module containing useful plotting abstractions on top of matplotlib. """

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec

class Plot2D:
    'Pseudocolor plot of a 2D array with optional 1D slices attached.'

    def __init__(self, arr2d, h_axis, v_axis, xlabel=r'', ylabel=r'', zlabel=r''):
        r"""
        >>> uu = np.linspace(0, np.pi, 128)
        >>> data = np.cos(uu - 0.5) * np.cos(uu.reshape(-1, 1) - 1.0)
        >>> p2d = Plot2D(data, uu, uu, xlabel=r'$x$ ($\mu$m)',
                ylabel=r'$y$ ($\mu$m)', zlabel=r'$\rho$ (cm$^{-3}$)')
        """
        self.data = arr2d
        self.h_axis = h_axis
        self.v_axis = v_axis
        self.label = {'x':xlabel, 'y':ylabel, 'z':zlabel}
        self.extent = [np.min(h_axis), np.max(h_axis),
                       np.min(v_axis), np.max(v_axis)]
        self.param = {'figsize':(6.4, 6.4), 'cmap':'coolwarm'}


    @staticmethod
    def colorbar(mappable):
        r"""Constructs a scaled colorbar for a given plot.
        
        Parameters
        ----------
        mappable : The Image, ContourSet, etc. to which the colorbar applies.
        """
        ax = mappable.axes
        fig = ax.figure
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        return fig.colorbar(mappable, cax=cax)


    def cbar(self, ax=None, **kwargs):
        r"""
        >>> p2d.cbar(figsize=(4.8, 4.8))
        """
        if not ax:
            self.param.update(kwargs)
            fig, ax = plt.subplots(figsize=self.param['figsize'])

        im = ax.imshow(self.data, origin='lower',
                             extent = self.extent,
                             aspect='auto',
                             interpolation='none',
                             cmap=self.param['cmap'])
    
        ax.set_xlabel(self.label['x'])
        ax.set_ylabel(self.label['y'])
    
        if not ax:
            cbar = self.colorbar(im)
            cbar.set_label(self.label['z'])    



    def slice(self, hslice_pos=None, vslice_pos=None, **kwargs):
        r"""
        >>> p2d.slice(figsize=(6.4, 6.4), hslice_pos=100, vslice_pos=10,
                         vslice_opts={'ls':'--', 'color':'green'},
                         hslice_opts={'ls':'-', 'color':'red'})
        """
        slice_opts = {'ls': '--', 'color': '0.25'}
        hslice_opts = kwargs.get('hslice_opts', slice_opts)
        vslice_opts = kwargs.get('vslice_opts', slice_opts)

        if not hslice_pos:
            hslice_pos = self.data.shape[0] // 2 
        if not vslice_pos:
            vslice_pos = self.data.shape[1] // 2 
        if not vslice_opts:
            vslice_opts = slice_opts
    
        # A lot of information about GridSpec can be found here:
        # https://matplotlib.org/tutorials/intermediate/gridspec.html
        
        # Figure with three axes of different width and height. Axis-sharing
        # is used to make easier interactive exploration but could perfectly
        # be removed if not needed.
    
        # set up a GridSpec grid and add the desired Axes one by one
        self.param.update(kwargs)

        fig = plt.figure(figsize=self.param['figsize'])
        gs = GridSpec(2, 2, height_ratios=[1, 3], width_ratios=[3, 1])
        ax0 = fig.add_subplot(gs[1, 0])
        axh = fig.add_subplot(gs[0, 0], sharex=ax0)
        axv = fig.add_subplot(gs[1, 1], sharey=ax0)
        
        # Imshow global view. Some documentation about imshow:
        # https://matplotlib.org/gallery/images_contours_and_fields/image_demo.html
        self.cbar(ax0)
        
        # Just to show where the slices are done. Another example in Matplotlib
        # gallery about the relevant functions:
        # https://matplotlib.org/gallery/subplots_axes_and_figures/axhspan_demo.html
        ax0.axhline(y=self.v_axis[hslice_pos], **hslice_opts)  ##----##
        ax0.axvline(x=self.h_axis[vslice_pos], **vslice_opts)  ## | ##
        # One can also annotate the slice markers. More information
        # about `ax.annotate` here:
        # https://matplotlib.org/tutorials/text/annotations.html
        ax0.annotate('{:07.3f}'.format(self.v_axis[hslice_pos]), xy=(self.h_axis[3], self.v_axis[hslice_pos+3]), xycoords='data', color=hslice_opts['color'])
        ax0.annotate('{:07.3f}'.format(self.h_axis[vslice_pos]), xy=(self.h_axis[vslice_pos-5], self.v_axis[self.v_axis.shape[0]//2 - 5]), xycoords='data', color=vslice_opts['color'], rotation='vertical')
        
        # Horizontal Slice Profile
        axh.set_xmargin(0)  # otherwise ax0 may have white margins
        axh.set_ylabel(self.label['z'])
        axh.plot(self.h_axis, self.data[hslice_pos], **hslice_opts)
        #axh.set_ylim(-1, 1)
        #axh.set_yticks([-1, 0, 1])
        # NB: More examples of fancy customisation of ticks and formatters here:
        # https://matplotlib.org/gallery/ticks_and_spines/tick-locators.html
        # https://matplotlib.org/gallery/ticks_and_spines/tick-formatters.html
        
        # Vertical Slice Profile
        axv.set_ymargin(0)  # otherwise ax0 may have white margins
        axv.set_xlabel(self.label['z'])
        axv.plot(self.data[vslice_pos], self.v_axis, **vslice_opts)
        #axv.set_xlim(-1, 1)
        #axv.set_xticks([-1, 0, 1])
        
        # "Despine" the slice profiles and hide the relevant axis
        for ax, spines in ((axh, ('top', 'bottom', 'right')),
                           (axv, ('top', 'left', 'right'))):
            for sp in spines:
                ax.spines[sp].set_visible(False)
        axh.xaxis.set_visible(False)
        axv.yaxis.set_visible(False)
        
        # Tweak a bit the figure layout
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.03, hspace=0.03)
    
        #return fig
        
    


def plot1d_break_x(fig, h_axis, v_axis, param, slice_opts):
    r"""
    >>> uu = np.linspace(0, np.pi, 128)
    >>> data = np.cos(uu - 0.5) * np.cos(uu.reshape(-1, 1) - 1.0)
    >>> fig, ax = plt.subplots(figsize=(8, 3.2))
    >>> plot1d_break_x(fig, uu, data[data.shape[0]//2], {'xlim_left':(0,1), 'xlim_right':(2,3),
        'xlabel':r'$x$ ($\mu$m)', 'ylabel':r'$\rho$ (cm$^{-3}$)'}, {'ls': '--', 'color': '0.5'})
    >>> fig
    """
    ax_left = fig.axes[0]
    divider = make_axes_locatable(ax_left)
    ax_right = divider.new_horizontal(size="100%", pad=1)
    fig.add_axes(ax_right)

    ax_left.plot(h_axis, v_axis, **slice_opts)
    ax_left.set_ylabel(param['ylabel'])
    ax_left.set_xlabel(param['xlabel'])

    ax_left.set_xlim(*param['xlim_left'])
    ax_left.yaxis.tick_left()
    ax_left.tick_params(labelright='off')
    ax_left.spines['right'].set_visible(False)

    ax_right.plot(h_axis, v_axis, **slice_opts)
    ax_right.set_ylabel(param['ylabel'])
    ax_right.set_xlabel(param['xlabel'])
    ax_right.yaxis.set_label_position("right")
    
    ax_right.set_xlim(*param['xlim_right'])
    ax_right.yaxis.tick_right()
    ax_right.spines['left'].set_visible(False)

    # From https://matplotlib.org/examples/pylab_examples/broken_axis.html
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax_left.transAxes, color='k', clip_on=False)
    ax_left.plot((1-d,1+d), (-d,+d), **kwargs)
    ax_left.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax_right.transAxes)  # switch to the right axes
    ax_right.plot((-d,+d), (1-d,1+d), **kwargs)
    ax_right.plot((-d,+d), (-d,+d), **kwargs)
