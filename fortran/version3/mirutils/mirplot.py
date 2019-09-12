# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
"""

author
======
  Novimir Antoniuk Pablant
    - npablant@pppl.gov
    - novimir.pablant@amicitas.com


description
===========
An interface to matplotlib that allows specifiction of complex plots
though parameter dictionaries.


example
=======

There are two ways to use mirplot:
  - quickplot function
  - MirPlot object

The real power of mirplot is in the MirPlot object, however the quickplot 
function is convenient for quickly creating simple plots.


quickplot
---------

The input for quickplot is a dictionary, or a list of dictionaries. 
Each dictionary contains a plot to add.

The simplest example:

  import numpy as np
  import mirplot

  x = np.arange(10)
  y = x
  mirplot.quickplot({'x':x, 'y':y})

Any supported plot properties can be added to the plot dictionary:

  mirplot.quickplot({'x':x, 'y':y, 'ybound':[0.0, 10.0], 'color':'red'})

To plot to a pdf file add a few keywords:

  mirplot.quickplot({'x':x, 'y':y}, pdf=True, filename='temp.pdf')

To add multiple subplots to a single figure we simple add a list of plot 
dictionaries:

  mirplot.quickplot([{'x':x, 'y':y}, {'x':x, 'y':y**2}])

To add multiple lines to a single subplot we need to include a plot name.

  mirplot.quickplot([{'name':'plot one', 'x':x, 'y':y}, {'name':'plot one', 'x':x, 'y':y**2}])

Each unique (or empty) plot name will results in a new subplot in the figure.


MirPlot
-------

The mirplot object is a list type where each element is a plot dictionary.

The simplest example:

  import numpy as np
  import mirplot

  x = np.arange(10)
  y = x

  plotobj = mirplot.MirPlot()
  plotobj.append({'x':x, 'y':y})
  plotobj.plotToScreen()




supported plot properties
=========================

type (string)
  Allowed Values: line, errorbar, scatter, fill_between, hline, vline, hspan, vspan

xerr (array)
yerr (array)

xbound (string)
ybound (string)

xlabel (string)
ylabel (string)

label (string)
  A label to use in the legend, if the the legend is active.

legend (bool)
  Set to true to show the legend in this subplot.

color (tuple or string)
  The color for the given plot.

colorscheme (string)
  Chose a color scheme for automatic color generation.
colorscheme_index (string)
  The index of the color scheme to use for the current plot.

alpha
markerfacecolor
markeredgecolor
facecolor

supress_xticklabels (bool)
  If true the xticklabels will be suppressed.

linestyle
  A string specifing the line style.

linewidth
  A value specifiing the line thickness.

marker
  A string specifing the marker style.

aspect
  A string or number specfiing the aspect ratio.
  
spported figure properties
==========================

figure_title (string)
figure_xlabel (string)

single_xlabel (bool)
  Only place an xlabel on the bottom most plot.  Suppress all other xlabels.
  The label will be taken from 'figure_xlabel' if given, or from the last 
  plot otherwise.

single_xticklabels (bool)
  Only place an xticklabel on the bottommost plot.  Suppress all other xticklabels.
  The tickmarks will be taken the last plot.

"""

import logging
import matplotlib
import numpy as np
from collections import OrderedDict

import mircolors

def quickplot(plotinput=None, filename=None, **keywords):
    """Quickly make a simple plot from a dictionary."""

    if plotinput is None:
        plotinput = [keywords]
    elif isinstance(plotinput, dict):
        plotinput = [plotinput]
                  
    mirplot = MirPlot()

    for plotdict in plotinput:
        mirplot.append(plotdict)

    if not keywords.get('pdf', False):
        mirplot.plotToScreen()
    else:
        mirplot.plotToPdf(filename)

    return mirplot


def getSchemeColor(name, index_input, num_colors=None):
    """
    Return a color for the current scheme.

    return
    ------

    A 4 element tuple is returned:
    (red, green, blue, alpha)


    inputs
    ------

    color_scheme (string)
      A name for the color scheme to use.

    index (int)
      The index of the color to take from the color scheme.

    keywords
    --------

    num_colors (int)
      The number of colors needed as part of the scheme. 
      This only affect gradient based schemes.

    """
    
    if num_colors is None:
        num_colors = 5

    # Make the index cycle if it exceeds the number of colors.
    index = index_input % 5

    gradient = mircolors.getColorGradient(
        mircolors.Normalize(0.0, float(num_colors-1))
        ,name)

    return gradient(float(index))


class MirPlot(list):
    """
    An interface to matplotlib that allows specifiction of complex plots
    though parameter dictionaries.

    See the module doc string for detailed documentation.
    """

    all_figure_properties = [
        'figure_title'
        ,'single_xticklabels'
        ,'figure_xlabel'
        ,'single_xlabel']
    
    all_subfig_properties = [
        'xerr'
        ,'yerr'
        ,'title'
        ,'xbound'
        ,'ybound'
        ,'xlim'
        ,'ylim'
        ,'xlabel'
        ,'ylabel'
        ,'xscale'
        ,'yscale'
        ,'tick_params'
        ,'supress_xticklabels'
        ,'legend'
        ,'legend_location'
        ,'legend_fontsize'
        ,'aspect']

    all_plot_properties = [
        'linestyle'
        ,'linewidth'
        ,'marker'
        ,'color'
        ,'colorscheme'
        ,'alpha'
        ,'label'
        ,'markerfacecolor'
        ,'markeredgecolor'
        ,'markersize'
        ,'facecolor']
    
    def __init__(self, *args, **kwargs):
        super(MirPlot, self).__init__(*args, **kwargs)

        self.log = logging.getLogger(self.__class__.__name__)
                
        self.user_properties = {}
        self.properties = {}

        self.subfig_properties = {}


    def plotToScreen(self):
        """
        Plot to screen using the default backend.

        This works best if run from within ipython with pylab.
        """

        from matplotlib import pyplot

        self.figure = pyplot.figure(figsize=self._getFigureSize())
        self.addPlotsToFigure()

        if not 'inline' in matplotlib.get_backend():
            self.figure.show()

        
    def plotToPdf(self, filename=None):
        """Plot to PDF."""

        if filename is None:
            raise NotImplementedError('A output filename is currently required.')

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_pdf import FigureCanvasPdf
        
        self.figure = Figure(figsize=self._getFigureSize())
        self.canvas = FigureCanvasPdf(self.figure)

        self.addPlotsToFigure()

        self.figure.savefig(filename)

        self.log.info('Saved figure to file: {}'.format(filename))


    def _autonamePlots(self):
        """
        Automatically name any plots that were not given a name by the user.
        """

        count = 0
        for ff in self:
            if not 'name' in ff:
                ff['name'] = '_autoname_{:02d}'.format(count)
                count += 1


    def _getSubfigNames(self, autoname=False):
        """
        Get a list of the subfigure names.
        """

        namelist = []

        for ff in self:
            if not 'name' in ff:
                if autoname:
                    # Autoname has not yet been called.  Call it now.
                    self._autonamePlots()
                else:
                    raise Exception('Some plots are not named.  Set autoname=True or call  _autonamePlots().')

            if not ff['name'] in namelist:
                namelist.append(ff['name'])

        return namelist

                
    def _countNumSubfigs(self):
        """
        Count the number of subfigures.
        """

        namelist = self._getSubfigNames(autoname=True)
        return len(namelist)


    def _getFigureSize(self):
        """
        Return the default figure size.

        The standard width is 6 inches.  The standard height is 2 inches for every subplot.

        return
        ------
          (width, height)
              The figure size in inches.
        """

        num_subfigs = self._countNumSubfigs()

        figure_width = 8
        figure_height = max(6, min(num_subfigs*3, 10))
        
        return (figure_width, figure_height)
        
        
        
    def addPlotsToFigure(self):
        """Add the plots to the figure."""

        # Go through the plot properties and cleaup if necessary.
        for ff in self:
            self.cleanupPlotProperties(ff)
            
        # Generate names if needed.
        self._autonamePlots()

            
        # Set default values for all of the plots.
        for ff in self:
            self.setPlotDefaults(ff)

        
        # Create the necessary subfigures. 
        self.subfigs = OrderedDict()
        for ff in self:
            if not ff['name'] in self.subfigs:
                self.subfigs[ff['name']] = {'plots':[]
                                            ,'properties':{}}

        num_subfigs = len(self.subfigs)

        # Add each plot to the appropriate subfigure.
        for ii, (key, subfig) in enumerate(self.subfigs.items()):
            # Create a new subplot
            subfig['ax'] = self.figure.add_subplot(num_subfigs, 1, ii+1)
            
            for ff in self:
                if ff['name'] == key:
                    subfig['plots'].append(ff)


        # Extract any properties relating to the overal figure.
        self.extractFigurePropertiesPlots()
        self.extractFigurePropertiesUser()
         
        for name, subfig in self.subfigs.items():
            self.addPlotsToSubfig(subfig)
            self.extractSubfigPropertiesPlots(subfig)
            self.extractSubfigPropertiesUser(name)
            self.applySubfigProperties(subfig)

        # Apply any properties relating to the overal figure.
        # This will overide any individual settings given
        # in the plots.
        self.applyFigureProperties()

        
    def cleanupPlotProperties(self, plot):
        """
        Check the plot properties and cleanup or provides errors.
        """
        if not 'x' in plot: plot['x'] = None
        if not 'y' in plot: plot['y'] = None
        if not 'xerr' in plot: plot['xerr'] = None
        if not 'yerr' in plot: plot['yerr'] = None

        # The capsize default is currently broken in matplotlib v2.0
        # explicityly set a default size here.
        if not 'capsize' in plot: plot['capsize'] = 2.0

        if 'x' in plot and plot['x'] is not None:
            if np.isscalar(plot['x']):
                plot['x'] = np.array([plot['x']])
            else:
                plot['x'] = np.array(plot['x'])
                
        if 'y' in plot and plot['y'] is not None:
            if np.isscalar(plot['y']):
                plot['y'] = np.array([plot['y']])
            else:
                plot['y'] = np.array(plot['y'])

        
    def setPlotDefaults(self, plot):
        plot.setdefault('type', 'line')

        if plot['x'] is None:
            plot['x'] =  np.arange(len(plot['y']))

        plot.setdefault('s', 15)

        
    def addPlotsToSubfig(self, subfig):
        """
        For a given subfigure.  Loop the though the plot dictionaries
        and add them to the subfigure axis.
        """
        
        axis = subfig['ax']

        # Create a dictionary to keep track of the current index for each colorscheme.
        # The idea here is that we track the index for each colorscheme separately
        # so that we can plot different data with differnt schemes on a single plot.
        colorscheme_indexes = {}
        colorscheme_num_colors = {}
        for plot in subfig['plots']:
            if 'colorscheme' in plot:
                if not plot['colorscheme'] in colorscheme_indexes:
                    colorscheme_indexes[plot['colorscheme']] = 0
                    colorscheme_num_colors[plot['colorscheme']] = 0
                    
                colorscheme_num_colors[plot['colorscheme']] += 1

                
        for plot in subfig['plots']:
            if plot['type'] == 'line':
                plotobj, = axis.plot(plot['x'], plot['y'])
            elif plot['type'] == 'errorbar':
                plotobj = axis.errorbar(
                    plot['x']
                    ,plot['y']
                    ,xerr=plot['xerr']
                    ,yerr=plot['yerr']
                    ,fmt='none'
                    ,capsize=plot['capsize'])
            elif plot['type'] == 'scatter':
                plotobj = axis.scatter(plot['x'], plot['y'], s=plot['s'], marker=plot.get('marker', None))
            elif plot['type'] == 'fill_between' or plot['type'] == 'fillbetween' :
                plotobj = axis.fill_between(plot['x'], plot['y'], plot['y1'])
            elif plot['type'] == 'hline':
                plotobj = axis.axhline(plot['y'][0])
            elif plot['type'] == 'vline':
                plotobj = axis.axvline(plot['x'][0])
            elif plot['type'] == 'hspan':
                plotobj = axis.axhspan(plot['y'][0], plot['y'][1])
            elif plot['type'] == 'vspan':
                plotobj = axis.axvspan(plot['x'][0], plot['x'][1])
            else:
                raise Exception('Plot type unknown: {}'.format(plot['type']))

            if 'color' in plot:
                color = plot['color']
            elif 'colorscheme' in plot:
                # Retrieve the next color in the color scheme, and update the
                # colorscheme counter.
                color = getSchemeColor(plot['colorscheme']
                                       ,colorscheme_indexes[plot['colorscheme']]
                                       ,colorscheme_num_colors[plot['colorscheme']])
                colorscheme_indexes[plot['colorscheme']] += 1
            else:
                color = None

            # Certain plot types actually consist of collections of line objects.
            # Loop through all objects and set the appropriate properties.
            if plot['type'] == 'errorbar':
                line_list = []
                line_list.extend(plotobj[1])
                line_list.extend(plotobj[2])

            else:
                line_list = [plotobj]

            for line in line_list:
                if 'linestyle' in plot:
                    line.set_linestyle(plot['linestyle'])
                    
                if 'linewidth' in plot:
                    line.set_linewidth(plot['linewidth'])
                    
                if 'marker' in plot:
                    if plot['type'] == 'scatter':
                        pass
                    else:
                        line.set_marker(plot['marker'])

                if color is not None:
                    line.set_color(color)

                if 'alpha' in plot:
                    line.set_alpha(plot['alpha'])

                if 'label' in plot:
                    line.set_label(plot['label'])

                if 'markerfacecolor' in plot:
                    line.set_facecolor(plot['markerfacecolor'])

                if 'markeredgecolor' in plot:
                    line.set_edgecolor(plot['markeredgecolor'])

                if 'markersize' in plot:
   
                    if plot['type'] == 'scatter':
                        sizes = line.get_sizes()
                        sizes[:] = plot['markersize']
                        line.set_sizes(sizes)
                    else:
                        line.set_markersize(plot['markersize'])

                if 'facecolor' in plot:
                    line.set_facecolor(plot['facecolor'])
                
            
    def applySubfigProperties(self, subfig):

        # Choose default plot properties if needed.
        if not 'ybound' in subfig['properties']:
            yrange = np.array([np.nanmin(subfig['plots'][0]['y']), np.nanmax(subfig['plots'][0]['y'])])
            # Note: I am doing this twice now for the first plot.
            for plot in subfig['plots']:
                yrange[0] = np.nanmin([yrange[0], np.nanmin(plot['y'])])
                yrange[1] = np.nanmax([yrange[1], np.nanmax(plot['y'])])

            subfig['properties']['ybound'] = yrange + np.array([-0.1, 0.1])*(yrange[1]-yrange[0])
        if not 'legend_fontsize' in subfig['properties']:
            subfig['properties']['legend_fontsize'] = 12.0
        if not 'legend_framealpha' in subfig['properties']:
            subfig['properties']['legend_framealpha'] = 0.7

        
        # Set the axis properties.
        # In some cases the order of these statements is important.
        # (For example xscale needs to come before xbound I think.)
        axis = subfig['ax']
        properties = subfig['properties']

        axis.tick_params('both'
                        ,direction='in'
                        ,which='both'
                        ,top='on'
                        ,bottom='on'
                        ,left='on'
                        ,right='on')
        
        if 'xscale' in properties:
            if properties['xscale'] == 'log':
                nonposx = 'clip'
            else:
                nonposx = None
            axis.set_xscale(properties['xscale'], nonposx=nonposx)
        
        if 'yscale' in properties:
            if properties['yscale'] == 'log':
                nonposy = 'clip'
            else:
                nonposy = None
            axis.set_yscale(properties['yscale'], nonposy=nonposy)
        if 'title' in properties:
            axis.set_title(properties['title'])
        if 'xbound' in properties:
            axis.set_xbound(properties['xbound'])
        if 'ybound' in properties:
            axis.set_ybound(properties['ybound'])
        if 'xlim' in properties:
            axis.set_xlim(properties['xlim'])
        if 'ylim' in properties:
            axis.set_ylim(properties['ylim'])
        if 'xlabel' in properties:
            axis.set_xlabel(properties['xlabel'])
        if 'ylabel' in properties:
            axis.set_ylabel(properties['ylabel'])
        if 'supress_xticklabels' in properties:
            axis.set_xticklabels([])
        if 'legend' in properties:
            if properties['legend']:
                axis.legend(loc=properties.get('legend_location')
                            ,fontsize=properties.get('legend_fontsize')
                            ,framealpha=properties.get('legend_framealpha')
                           )
        if 'aspect' in properties:
            axis.set_aspect(properties['aspect'])
                
    def extractSubfigPropertiesPlots(self, subfig):

        # First see if the user defined any of the plot properties.
        for plot in subfig['plots']:
            for key in self.all_subfig_properties:
                if key in plot:
                    subfig['properties'][key] = plot[key]

                
    def extractSubfigPropertiesUser(self, name):
        subfig = self.subfigs[name]

        if name in self.subfig_properties:
            for key in self.subfig_properties[name]:
                if key in self.all_subfig_properties:
                    subfig['properties'][key] = self.subfig_properties[name][key]
                else:
                    raise Exception('Subfigure property not valid: {}'.format(key))

                
    def setPlotProperties(self, name, **properties):
        if not name in self.subfig_properties:
            self.subfig_properties[name] = {}
            
        for key in properties:
                self.subfig_properties[name][key] = properties[key]      
                

    def extractFigurePropertiesPlots(self):
        for plot in self:
            for prop in self.all_figure_properties:
                if prop in plot:
                    self.properties[prop] = plot[prop]


    def extractFigurePropertiesUser(self):
        for key in self.user_properties:
            if key in self.all_figure_properties:
                self.properties[key] = self.user_properties[key]
            else:
                raise Exception('Figure property not valid: {}'.format(key))

                
    def applyFigureProperties(self):
        
        # NOTE: The order of parsing here is important as some properties modify
        #       other properties.

        if 'figure_title' in self.properties:
            self.figure.suptitle(self.properties['figure_title'])

        if 'single_xticklabels' in self.properties:
            if self.properties['single_xticklabels']:
                keys_list = list(self.subfigs.keys())
                # Suppres the x tick labels for every subfigure except the last one.
                for ii in range(len(self.subfigs)-1):
                    self.subfigs[keys_list[ii]]['ax'].set_xticklabels([])

        if 'figure_xlabel' in self.properties:
            # Allows a figure wide xlabel to be applies to the plot.
            # This will suppress all individual plot labels.
            keys_list = list(self.subfigs.keys())
            self.subfigs[keys_list[-1]]['ax'].set_xlabel(self.properties['figure_xlabel'])

            self.properties['single_xlabel'] = True

        if 'single_xlabel' in self.properties:
            if self.properties['single_xlabel']:
                keys_list = list(self.subfigs.keys())
                # Suppres the x labels for every subfigure except the last one.
                for ii in range(len(self.subfigs)-1):
                    self.subfigs[keys_list[ii]]['ax'].set_xlabel('')


    def setFigureProperties(self, *args, **kwargs):
        if args:
            properties = args[0]
        else:
            properties = kwargs

        for key in properties:
                self.user_properties[key] = properties[key]          

                
    def getPlot(self, name):
        return self.subfigs[name]
