## Ideal Guiding Model

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import scipy.linalg

class Guider(object):
    """Calculates optimum guider corrections.
    """
    def __init__(self, x0, y0, x, y):
        """
        Find the scale, rotation and offsets that minimizes the residuals between x,y and x0,y0.
        """
        assert x0.shape == y0.shape, 'x0,y0 have different shapes.'
        assert x.shape == y.shape, 'x,y have different shapes.'
        assert len(x.shape) == 2, 'x,y have unexpected shape.'
        assert x.shape[0] == x0.shape[0], 'x,y and x0,y0 have different lengths.'

        nxy, nt = x.shape
        A = np.empty((2 * nxy, 4))
        self.scale = np.empty((nt,))
        self.rotation = np.empty((nt,)) * u.rad
        self.dx = np.empty((nt,)) * u.m
        self.dy = np.empty((nt,)) * u.m
        self.nt = nt
        self.platescale = 217.7358 * u.mm / u.deg
        
        xy0 = np.concatenate([x0.si.value.flat, y0.si.value.flat])
        for it in range(nt):
            xt = x[:, it].si.value.flatten()
            yt = y[:, it].si.value.flatten()
            xyt = np.concatenate([xt, yt])
            for i in range(nxy):
                A[i, :] = (xt[i], -yt[i], 1. ,0.)
                A[nxy + i, :] = (yt[i], xt[i], 0., 1.)
            params, xy_residuals, rank, sing = scipy.linalg.lstsq(A, xy0)
            
            self.scale[it] = np.sqrt(params[0]**2 + params[1]**2)
            self.rotation[it] = np.arctan2(params[1], params[0]) * u.rad
            self.dx[it] = params[2] * u.m
            self.dy[it] = params[3] * u.m
            
        self.scos = self.scale * np.cos(self.rotation)
        self.ssin = self.scale * np.sin(self.rotation)
        
        self.x0 = x0
        self.y0 = y0
        self.x_before = x
        self.y_before = y
        self.x_after, self.y_after = self.correct(x, y)

    def correct(self, x, y):
        """Apply guiding corrections to the focal plane positions x, y.
        """
        assert x.shape == y.shape, 'x,y have different shapes.'
        assert x.shape[-1] == self.nt, 'x,y have unexpected shape.'
        
        return self.scos * x - self.ssin * y + self.dx, self.ssin * x + self.scos * y + self.dy
    
    def plot(self, tai, zoom=3000, field_radius=None, fiber_radius=None, save=None):
        """
        """
        plt.figure(figsize=(12, 12))
        assert len(tai.shape) == 1 and len(tai) == self.nt, 'tai has unexpected shape.'
        
        assert field_radius is not None
        rmax = field_radius.to(u.mm).value            
        plt.xlim(-1.01 * rmax, +1.01 * rmax)
        plt.ylim(-1.01 * rmax, +1.01 * rmax)
        plt.axis('off')
        outline = plt.Circle((0,0), rmax, edgecolor='black', facecolor='none')
        plt.gca().add_artist(outline)
        plt.gca().set_aspect(1.0)
        
        # Draw the nominal guide hole centers.
        plt.scatter(self.x0.to(u.mm).value, self.y0.to(u.mm).value, marker='+', s=100, color='k')

        for i, (x0, y0) in enumerate(zip(self.x0, self.y0)):
            if fiber_radius is not None:
                fiber = plt.Circle(
                    (x0.to(u.mm).value, y0.to(u.mm).value),
                    zoom * fiber_radius.to(u.mm).value, edgecolor='k', facecolor='none', ls='dotted')
                plt.gca().add_artist(fiber)
            plt.plot((x0 + zoom * (self.x_before[i] - x0)).to(u.mm).value,
                     (y0 + zoom * (self.y_before[i] - y0)).to(u.mm).value, 'b-',
                     label=(None if i else 'No Guiding'))
            plt.plot((x0 + zoom * (self.x_after[i] - x0)).to(u.mm).value,
                     (y0 + zoom * (self.y_after[i] - y0)).to(u.mm).value, 'r-',
                    label=(None if i else 'Ideal Guiding'))
        plt.legend(loc='lower left')
        
        plt.tight_layout()
        if save:
            plt.savefig(save)

if __name__ == '__main__':
    #guider = Guider(guide_x0, guide_y0, guide_x, guide_y)
    pass