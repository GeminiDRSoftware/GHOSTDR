from __future__ import division, print_function

import sys
import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Plotter:
    def __init__(self, fname1, fname2):
        """ A class to handle the plotting and animation of two spectra.
            The first is the extracted spectrum produced by the pipeline,
            and the second is the check spectrum produced by the simulator.
        """

        self.fname1 = fname1
        self.fname2 = fname2

        # Open the science frame
        scihdu = pf.open(self.fname1)
        scidata = scihdu[1].data

        # Gumby check for hires
        ishigh = False
        if 'high' in self.fname1:
            ishigh = True

        # Open the check frame
        chkhdu = pf.open(self.fname2)
        chkdata = chkhdu[0].data

        self.fig = plt.figure()

        # Accumulate the data
        # This code evolved via a lot of experiments with different
        # styles of plots and could probably be made a lot tidier but
        # I'm not spending time on it now.
        verts1 = []
        verts2 = []
        xs = np.arange(len(scidata[0,:,0]))
        self.zs = np.arange(scidata.shape[0])
        for i in self.zs:
            chkorder = chkdata[i]
            sciorder = scidata[i,:,0]
            # Each order of the check spectrum is scaled so that
            # its median matches the median of the extracted spectrum.
            # This is because the check spectrum is in arbitrary units.
            chkmedian = np.median(chkorder)
            if abs(chkmedian) > 0.01:
                scale = np.median(sciorder) / np.median(chkorder)
            else:
                scale = 1
            verts1.append(list(zip(xs, sciorder)))
            verts2.append(list(zip(xs, scale * chkorder)))

        self.v1 = np.asarray(verts1)
        self.v2 = np.asarray(verts2)

        # Work out the ranges for the axes
        x0 = self.v1[:,:,0].min()
        x1 = self.v1[:,:,0].max()
        # Scale the y axis so we screen out the most extreme points
        y0, y1 = np.percentile(scidata[:,:,0], (2,98))

        # Construct the two (empty) plots
        plt.subplot(2, 1, 1)
        plt.title('%s' % self.fname1)
        self.l1, = plt.plot([], [])
        plt.xlim(x0, x1)
        plt.ylim(y0, y1)
        plt.grid()

        plt.subplot(2, 1, 2)
        plt.title('%s' % self.fname2)
        self.l2, = plt.plot([], [])
        plt.xlim(x0, x1)
        plt.ylim(y0, y1)
        plt.grid()

        # Don't waste space, use a tight layout
        self.fig.tight_layout()

    def update(self, num):
        """ The function that allows us to animate the orders (i.e. extensions).
            It also lets us plot a single specific extension if we want to.
        """
        plt.subplot(2, 1, 1)
        plt.title('%s %d' % (self.fname1, num))
        self.l1.set_data(self.v1[num,:,0], self.v1[num,:,1])
        plt.subplot(2, 1, 2)
        plt.title('%s %d' % (self.fname2, num))
        self.l2.set_data(self.v2[num,:,0], self.v2[num,:,1])
        return self.l1, self.l2

    def onClick(self, event):
        """ Pause/Resume the animation, based on a click in the plot.
        """
        if self.anim_running:
            self.animator.event_source.stop()
            self.anim_running = False
        else:
            self.animator.event_source.start()
            self.anim_running = True

    def start_animation(self):
        """ Start the animation running.
        """
        # Setting blit=True for the animation just makes a mess on my screen
        self.animator = animation.FuncAnimation(self.fig, self.update, self.zs, fargs=None, interval=500)
        # The user can pause the animation by clicking on the plot
        self.anim_running = True
        self.fig.canvas.mpl_connect('button_press_event', self.onClick)
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) == 3:
        # Animate all extensions
        plotter = Plotter(sys.argv[1], sys.argv[2])
        plotter.start_animation()
    elif len(sys.argv) == 4:
        # Just show the given extension
        plotter = Plotter(sys.argv[1], sys.argv[2])
        plotter.update(int(sys.argv[3]))
        plt.show()
    else:
        print("Usage: %s extracted.fits check.fits [ext_number]" % sys.argv[0])
