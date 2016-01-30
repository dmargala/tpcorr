#!/usr/bin/env python
"""
Fits tabulated throughput corrections to a model
"""
import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from sklearn import linear_model
import scipy.optimize

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", action="store_true",
        help="print verbose output")
    parser.add_argument("-o", "--output", type=str, default=None,
        help="output file base name")
    parser.add_argument("-i", "--input", type=str, default=None,
        help="required input file")
    parser.add_argument("--nexamples", type=int, default=-1,
        help="number of example fits to plot, (array slice 'end' value)")
    parser.add_argument("--seed", type=int, default=2810,
        help="random seed for example choices")
    parser.add_argument("--scipy", action="store_true",
        help="use scipy.optimize instead of sklearn")
    parser.add_argument("--exp-info", type=str, default=None,
        help="add exp info")
    parser.add_argument("--skip-plots", action="store_true",
        help="don't make plots")
    args = parser.parse_args()

    if args.verbose:
        print 'Reading file: %s' % args.input

    # for now, construct the title from the input filename
    title = '-'.join((args.input).split('-')[1:])[:-4]

    # The input data is text file where each line coresponds to 
    # a target's throughput correction vector
    data = np.loadtxt(args.input, ndmin=2)

    try:
        nentries, ntokens = data.shape
    except ValueError:
        print data.shape, args.input

    # the first two columns are xfocal and yfocal positions of the target
    nidtokens = 3
    # the rest are the tabulated throughput correction values
    npoints = ntokens - nidtokens

    # the throughput correction vectors span the range 3500A to 10500A
    xvalues = np.linspace(3500, 10500, npoints, endpoint=True)

    nparams = 3
    results = np.empty((nentries, nparams))
    chisqs = np.empty(nentries)

    # fit individual entries
    for i,row in enumerate(data):
        yvalues = row[nidtokens:]
        if args.scipy:
            # chisq function for our model
            def chisq(params):
                sigma = 1e-1
                x0 = np.exp(params[0])
                pred = 1+params[1]*np.log(xvalues/x0)+params[2]*np.log(xvalues/x0)**2
                residuals = (yvalues - pred)/sigma
                return np.dot(residuals,residuals)
            params0 = np.array([np.log(7000),1,-.5])
            result = scipy.optimize.minimize(chisq, params0, options={'maxiter':10000},
                method='Nelder-Mead')
            # save fit results
            results[i,:] = result.x
            results[i,0] = np.exp(result.x[0])
            chisqs[i] = result.fun
            if not result.success:
                # chisq function for our model
                def chisq(params):
                    sigma = 1e-1
                    x0 = params[0]
                    pred = 1+params[1]*np.log(xvalues/x0)+params[2]*np.log(xvalues/x0)**2
                    residuals = (yvalues - pred)/sigma
                    return np.dot(residuals,residuals)
                params0 = np.array([7000,1,-.5])
                result = scipy.optimize.minimize(chisq, params0, options={'maxiter':10000},
                    method='SLSQP', bounds=((1,None),(None,None),(None,None)))
                # save fit results
                results[i,:] = result.x
                chisqs[i] = result.fun
                if not result.success or result.x[0] == 1:
                    print 'failed on %s-%d: %s' % (title, row[0], results[i])
        else:
            # construct a matrix for the model A + B*Log(x) + C*(Log(x))^2
            xmatrix = np.ones((npoints, nparams))
            xmatrix[:,1] = np.log(xvalues)
            xmatrix[:,2] = np.log(xvalues)**2
            regr = linear_model.LinearRegression(fit_intercept=False)
            regr.fit(xmatrix, yvalues)
            a,b,c = regr.coef_
            # transform fit parameters to "physical" params
            x0 = np.exp(-0.5*(b-np.sqrt(b*b+4*(1-a)*c))/c)
            a1 = b + 2*c*np.log(x0)
            a2 = c
            # save fit results
            chisqs[i] = np.nansum((regr.predict(xmatrix)-yvalues)**2)
            results[i,:] = [x0, a1, a2]

    if args.verbose:
        print 'Mean fitted params (x0, a1, a2): (%.4g, %.4g, %.4g)' % tuple(np.nanmean(results, axis=0))
        print 'Mean chisq: %.4g' % np.mean(chisqs)
    else:
        print '%s %.4g %.4g %.4g %.4g' % tuple([title, np.mean(chisqs)] + np.nanmean(results, axis=0).tolist())

    # save results to file
    output = '%s-%s' % (args.output, '-'.join((args.input).split('-')[1:]))
    if args.verbose:
        print 'Saving results to: %s' % output 
    np.savetxt(output, results)

    if args.skip_plots:
        return 0

    if args.verbose:
        print 'Creating fit summary plot...'

    # save results summary plot
    fig = plt.figure(figsize=(12,12))

    ax1 = plt.subplot2grid((4,3), (0,0), colspan=3)
    ax2 = plt.subplot2grid((4,3), (1,0))
    ax3 = plt.subplot2grid((4,3), (2,1))
    ax4 = plt.subplot2grid((4,3), (3,2))

    ax5 = plt.subplot2grid((4,3), (2,0))
    ax6 = plt.subplot2grid((4,3), (3,0))
    ax7 = plt.subplot2grid((4,3), (3,1))

    ax8 = plt.subplot2grid((4,3), (1,2))

    # compare the raw throughput corrections with fit results
    plt.sca(ax1)
    # randomly select example fits to plot
    np.random.seed(args.seed)
    for params in np.random.permutation(results)[:args.nexamples]:
        yvalues = 1 + params[1]*np.log(xvalues/params[0]) + params[2]*(np.log(xvalues/params[0]))**2
        plt.plot(xvalues, yvalues, c='black', alpha=.2, lw=.05)
    # shade max-min region of the raw throughput correction vectors
    plt.fill_between(xvalues, data[:,nidtokens:].min(axis=0),data[:,nidtokens:].max(axis=0), alpha=.5, lw=0)
    # manually construct a legend for this plot
    blue_patch = mpatches.Patch(color='blue', alpha=.5, label=r'Prediction Coverage')
    black_line = mlines.Line2D([], [], color='black', alpha=.5, label=r'$1 + a_1 \log(x/x_0) + a_2 (\log(x/x_0))^2 $')
    plt.legend(handles=[blue_patch, black_line], loc=2)
    # add labels and set ranges
    plt.xlabel(r'Wavelength $(\AA)$')
    plt.ylabel('Throughput Correction')
    plt.ylim([0, 3])
    plt.xlim([3500, 10500])
    plt.title(title)
    plt.grid()        

    # plot the fit parameter distributions
    def plot_param_dist(params, binspec, color, xlabel):
        xmin, xmax, nxbins = binspec
        plt.hist(params, bins=np.linspace(xmin, xmax, nxbins+1), facecolor=color, alpha=.5)
        plt.ylabel('Counts')
        plt.xlabel(xlabel)
        plt.xlim([xmin, xmax])
        plt.grid()
        textstr = '$\mathrm{mean}=%.2f$\n$\mathrm{median}=%.2f$\n$\mathrm{std}=%.2f$' % (np.nanmean(params), np.nanmedian(params), np.nanstd(params))
        props = dict(boxstyle='round', facecolor='white')
        plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, va='top', ha='right', bbox=props)

    plt.sca(ax2)
    plot_param_dist(results[:,0], (3500,7000,50), 'blue', r'$x_0$')

    plt.sca(ax3)
    plot_param_dist(results[:,1], (-0.55,2.5,50), 'green', r'$a_1$')

    plt.sca(ax4)
    plot_param_dist(results[:,2], (-1,0.5,50), 'red', r'$a_2$')

    # plot the fit parameter distributions
    def plot_param_scatter(x, y, xlim, ylim, xlabel, ylabel):
        plt.plot(x, y, '+', ms=5)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.grid()
        # calculate correlation coefficient
        corr =  np.corrcoef(x,y)
        rho = corr[0,1]
        # add text box
        textstr = r'$\rho=%.2f$' % rho
        props = dict(boxstyle='round', facecolor='white')
        plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, va='top', ha='right', bbox=props)

    plt.sca(ax5)
    plot_param_scatter(results[:,0], results[:,1], (3500,7000), (-0.55,2.5), r'$x_0$', r'$a_1$')

    plt.sca(ax6)
    plot_param_scatter(results[:,0], results[:,2], (3500,7000), (-1,0.5), r'$x_0$', r'$a_2$')

    plt.sca(ax7)
    plot_param_scatter(results[:,1], results[:,2], (-0.55,2.5), (-1,0.5), r'$a_1$', r'$a_2$')

    plt.sca(ax8)
    ax8.axis('off')

    if args.exp_info:
        import json
        expinfo_filename = '%s-%s.json' % (args.exp_info, title)
        with open(expinfo_filename) as jsondata:
            expinfo = json.load(jsondata)
        keyfmt_pairs = [
            ('design_alt', '%.2f'),
            ('mean_alt', '%.2f'),
            ('design_ha', '%.2f'),
            ('mean_ha', '%.2f'),
            ('SEEING50', '%.2f'),
            ('mean_psf_fwhm', '%.2f')
        ]
        textstr = ''
        textstr += 'n_entries: %d\n' % nentries
        textstr += 'mean_chisq: %.4g\n' % np.mean(chisqs)
        textstr += 'nexp: %d\n' % len(expinfo['exposures'])
        textstr += '\n'.join([('%s: '+fmt) % (key, expinfo[key]) for key,fmt in keyfmt_pairs])
        plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, va='top', ha='right')

    plot_name = output[:-4]+'.png'
    if args.verbose:
        print 'Saving fit summary figure to file: %s' % plot_name
    fig.savefig(plot_name, bbox_inches='tight')

if __name__ == '__main__':
    main()
