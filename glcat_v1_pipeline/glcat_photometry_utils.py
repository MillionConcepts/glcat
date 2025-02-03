import numpy as np
import emcee
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Circle
from astropy.visualization import ZScaleInterval

# ncat = pd.read_csv('ncat_240616.csv',index_col=None)

def gaussian_flux_fraction(r, sigma):
    return 1 - np.exp(-r**2 / (2 * sigma**2))

def gaussian_flux_model(theta, r):
    total_flux, sigma, background = theta
    return total_flux * (1 - np.exp(-r**2 / (2 * sigma**2))) + background * np.pi * r**2

def log_likelihood(theta, r, flux, flux_err):
    model = gaussian_flux_model(theta, r)
    return -0.5 * np.sum(((flux - model) / flux_err)**2)

def log_prior(theta):
    total_flux, sigma, background = theta
    if (0 <= total_flux < 1000) and (0 < sigma < 10) and (0 <= background < 100):
        return 0.0
    return -np.inf
    
def log_probability(theta, r, flux, flux_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, r, flux, flux_err)

def mcmc_aperture_curve(aperture_radii,flux,flux_err,
                        nsteps = 1000, # number of MCMC steps
                        burnin = 200, # number of burn-in steps
                        nwalkers = 32, # number of MCMC walkers
                       ):
    ix = np.where(np.isfinite(flux))
    aperture_radii = aperture_radii[ix]
    flux = flux[ix]
    flux_err = flux_err[ix]

    # Set up the MCMC sampler
    ndim = 3  # Number of parameters to fit (total_flux, sigma, background)
    
    # Initialize the walkers
    #                total_flux,    sigma,   background
    initial_guess = [np.mean(flux), 5/2.355, (flux[-1]-flux[-2])/(np.pi * (aperture_radii[-1]**2-aperture_radii[-2]**2))]
    pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)
    
    # Run the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, 
                                    log_probability, 
                                    args=(aperture_radii, flux, flux_err),)
    burn_in = sampler.run_mcmc(pos, burnin)
    sampler.reset()
    sampler.run_mcmc(burn_in, nsteps)

    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    samples = sampler.get_chain()

    return flat_samples,samples

def get_percentile_ranges(flat_samples,percentiles=[16,50,84],labels = ["cps", "sigma", "bg_cps"],
                          quiet=True):
    assert len(percentiles) == 3
    if not quiet:
        print("Best-fit parameters:")
    stats = {l:[] for l in labels}
    for i in range(len(labels)):
        mcmc = np.percentile(flat_samples[:, i], percentiles)
        q = np.diff(mcmc)
        if not quiet:
            print(f"{labels[i]}: {mcmc[1]:.3f} +{q[1]:.3f} -{q[0]:.3f}")
        stats[labels[i]] = (mcmc[1],q[0],q[1])
    return stats

def plot_mcmc_results(flat_samples,aperture_radii,flux,flux_err,
                      percentiles=[16,50,84],
                      labels = ["cps", "sigma", "bg_cps"],
                      r = np.arange(1.5,17.3,0.1),
                      savefig=None,):

    # Compute the model predictions for each set of MCMC samples
    model_samples = np.zeros((len(flat_samples), len(r)))
    for i, theta in enumerate(flat_samples):
        model_samples[i] = gaussian_flux_model(theta, r)

    # Calculate the median and confidence intervals of the model predictions
    model_percentiles = np.percentile(model_samples, [16, 50, 84], axis=0)

    # Plot the data points
    plt.figure(figsize=(8, 6))
    plt.errorbar(aperture_radii, flux, yerr=flux_err, fmt='k.', label='data w/ 1-sigma error bars')

    # Plot the median model prediction and confidence intervals
    plt.plot(r, model_percentiles[1], 'b:', label='model w/ 68% conf. interval',alpha=0.5)
    plt.fill_between(r, model_percentiles[0], model_percentiles[2], color='b', alpha=0.2)#, label='68% confidence interval')
    #plt.semilogy()

    # Get percentile ranges and add to plot
    stats = get_percentile_ranges(flat_samples, percentiles, labels)
    # textstr = '\n'.join([f"{label}: {stats[label][0]:.3f} +{stats[label][2]:.3f} -{stats[label][1]:.3f}" for label in ['cps','bg_cps','sigma']])
    textstr = '\n'.join([f"{'flux'}: {stats['cps'][0]:.3f} +{stats['cps'][2]:.3f} -{stats['cps'][1]:.3f} cps     ",
                         f"{'bgnd'}: {stats['bg_cps'][0]:.3f} +{stats['bg_cps'][2]:.3f} -{stats['bg_cps'][1]:.3f} cps/as^2",
                         f"{'fwhm'}: {2.355*stats['sigma'][0]:.3f} +{2.355*stats['sigma'][2]:.3f} -{2.355*stats['sigma'][1]:.3f} as     "])

    plt.gcf().text(0.95, 0.15, textstr, fontsize=12, bbox=dict(facecolor='white', alpha=0.5), ha='right', va='bottom')

    plt.xlabel('aperture radius (as)')
    plt.ylabel('flux (cps)')
    plt.legend()
    plt.tight_layout()
    if savefig is not None:
        plt.savefig(savefig)
        plt.close()
    else:
        plt.show()

def plot_mcmc_walkers(samples,savefig=None):
    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    ndim = 3
    labels = ["cps", "sigma", "bg_cps"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    axes[-1].set_xlabel("step number");
    if savefig is not None:
        plt.savefig(savefig)
        plt.close()
    else:
        plt.show()

def plot_model_qa(samples,flux,flux_err,imgfn,
                    aperture_radii = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3]),
                      percentiles=[16,50,84],
                      labels = ["cps", "sigma", "bg_cps"],
                      r = np.arange(1.5,17.3,0.1),
                      savefig=None,):
    fig = plt.figure(figsize=(12,10))
    gs = GridSpec(6,6)

    # model plot
    mdl = fig.add_subplot(gs[:3,:3])
    mdl.errorbar(aperture_radii, flux, yerr=flux_err, fmt='k.', label='data w/ 1-sigma error bars')
    mdl.plot(r, model_percentiles[1], 'b:', label='model w/ 68% conf. interval',alpha=0.5)
    mdl.fill_between(r, model_percentiles[0], model_percentiles[2], color='b', alpha=0.2)#, label='68% confidence interval')
    stats = get_percentile_ranges(flat_samples, percentiles, labels)
    # textstr = '\n'.join([f"{label}: {stats[label][0]:.3f} +{stats[label][2]:.3f} -{stats[label][1]:.3f}" for label in ['cps','bg_cps','sigma']])
    textstr = '\n'.join([f"{'flux'}: {stats['cps'][0]:.3f} +{stats['cps'][2]:.3f} -{stats['cps'][1]:.3f} cps     ",
                        f"{'bgnd'}: {stats['bg_cps'][0]:.3f} +{stats['bg_cps'][2]:.3f} -{stats['bg_cps'][1]:.3f} cps/as^2",
                        f"{'fwhm'}: {2.355*stats['sigma'][0]:.3f} +{2.355*stats['sigma'][2]:.3f} -{2.355*stats['sigma'][1]:.3f} as     "])

    mdl.text(17.5, 0.01, textstr, fontsize=12, bbox=dict(facecolor='white', alpha=0.5), ha='right', va='bottom')

    mdl.set_xlabel('aperture radius (as)')
    mdl.set_ylabel('flux (cps)')
    mdl.legend()

    # walkers plot
    ndim = 3
    labels = ["cps", "sigma", "bg_cps"]
    for j in range(ndim):
        wkr = fig.add_subplot(gs[j,3:])
        wkr.plot(samples[:, :, j], "k", alpha=0.3)
        wkr.set_xlim(0, len(samples))
        wkr.set_title(labels[j])
        wkr.yaxis.set_label_coords(-0.1, 0.5)
        wkr.yaxis.tick_right()
        if j<2:
            wkr.set_xticks([])
    wkr.set_xlabel("step number");

    imgx = float(ncat.column('NUV_XCENTER')[source_ix].as_py())
    imgy = float(ncat.column('NUV_YCENTER')[source_ix].as_py())
    imsz = np.shape(ffull['CNT'])
    boxsz = 150
    x1, x2, y1, y2 = (max(int(imgy - boxsz), 0),
                    min(int(imgy + boxsz), imsz[0]),
                    max(int(imgx - boxsz), 0),
                    min(int(imgx + boxsz), imsz[1]))
    x1_,y1_=0,0

    ffull = pdr.read(imgfn)
    # full frame image
    ffi = fig.add_subplot(gs[3:,:3])
    ffi.imshow(ZScaleInterval()(ffull['CNT']/expt[0]),origin='lower',cmap='Greys')
    rect = Rectangle((y1 - y1_, x1 - x1_), 2 * boxsz, 2 * boxsz,
                    linewidth=1, edgecolor='y', facecolor='none',ls='solid')
    ffi.add_patch(rect)
    ffi.set_xticks([])
    ffi.set_yticks([])

    # sub frame image
    sfi = fig.add_subplot(gs[3:,3:])
    plt.imshow(ZScaleInterval()(ffull['CNT'][x1:x2,y1:y2]/expt[0]),origin='lower',cmap='Greys')
    circ1 = Circle((boxsz, boxsz), 17.5/1.5,
                linewidth=2, edgecolor='r', facecolor='none', ls='solid')
    sfi.add_patch(circ1)
    circ2 = Circle((boxsz, boxsz), 9.0/1.5,
                linewidth=1, edgecolor='r', facecolor='none', ls='dotted')
    sfi.add_patch(circ2)


    sfi.set_xticks([])
    sfi.set_yticks([])

    plt.tight_layout()
    if savefig is not None:
        plt.savefig(savefig)
        plt.close()
    else:
        plt.show()

# aperture_radii = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3])

# np.random.seed(1965)
# models = pd.DataFrame({})
# for i in tqdm(np.arange(len(ncat))):
#     flux = np.array([ncat[f'FUV_CPS_APER{a}'].iloc[i] for a,r in enumerate(aperture_radii)])
#     flux_err = np.array([ncat[f'FUV_CPS_ERR_APER{a}'].iloc[i] for a,r in enumerate(aperture_radii)])
#     ix = np.where(np.isfinite(flux))
#     flat_samples = fit_aperture_curve(aperture_radii[ix],flux[ix],flux_err[ix],nsteps=200)

#     a,b,c = (np.percentile(flat_samples[:, 0], [16, 50, 84])[1],
#              np.percentile(flat_samples[:, 1], [16, 50, 84])[1],
#              np.percentile(flat_samples[:, 2], [16, 50, 84])[1])
    
#     q_a = np.diff(np.percentile(flat_samples[:, 0], [16, 50, 84]))
#     q_b = np.diff(np.percentile(flat_samples[:, 1], [16, 50, 84]))
#     q_c = np.diff(np.percentile(flat_samples[:, 2], [16, 50, 84]))
    
#     models = pd.concat([models,
#                         pd.DataFrame({'cps':[a],
#                                       'cps_lower_conf':[q_a[0]],
#                                       'cps_upper_conf':[q_a[1]],
#                                       'sigma':[b],
#                                       'sigma_lower_conf':[q_b[0]],
#                                       'sigma_upper_conf':[q_b[1]],
#                                       'bg_cps':[c],
#                                       'bg_cps_lower_conf':[q_c[0]],
#                                       'bg_cps_upper_conf':[q_c[1]],
#                                       'bg_cps_annulus':[(flux[-1]-flux[-2])/(np.pi * (aperture_radii[-1]**2 - aperture_radii[-2]**2))],
#                                       })],
#                        axis=0)

# #    plt.figure()
# #    plt.errorbar(aperture_radii,flux,yerr=flux_err,fmt='k-',label='data')
# #    plt.plot(aperture_radii,
# #             gaussian_flux_model([a,b,c,], aperture_radii),'b:',
# #             label=f'flux,sigma,bg= ({a:.3f} {b:.3f} {c:.3f})')
# #    plt.legend()
# #    plt.semilogy();

# models.to_csv('ncat_fuv_fit_240616.csv',index=False)