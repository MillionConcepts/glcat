import numpy as np
import emcee
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

ncat = pd.read_csv('ncat_240616.csv',index_col=None)

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
    if (0 < total_flux < 1000) and (0 < sigma) and (0 < background < 100):
        return 0.0
    return -np.inf
    
def log_probability(theta, r, flux, flux_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, r, flux, flux_err)

def mcmc_aperture_curve(aperture_radii,flux,flux_err,
                        nsteps = 500, # number of MCMC steps
                        burnin = 200, # number of burn-in steps
                       ):
    ix = np.where(np.isfinite(flux))
    aperture_radii = aperture_radii[ix]
    flux = flux[ix]
    flux_err = flux_err[ix]

    # Set up the MCMC sampler
    ndim = 3  # Number of parameters to fit (total_flux, sigma, background)
    nwalkers = 32  # Number of MCMC walkers
    
    # Initialize the walkers
    #                total_flux,    sigma,   background
    initial_guess = [np.mean(flux), 5/2.355, (flux[-1]-flux[-2])/(np.pi * (aperture_radii[-1]**2-aperture_radii[-2]**2))]
    pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)
    
    # Run the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(aperture_radii, flux, flux_err))
    burn_in = sampler.run_mcmc(pos, burnin)
    sampler.reset()
    sampler.run_mcmc(burn_in, nsteps)

    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

    return flat_samples

def get_percentile_ranges(flat_samples,percentiles=[16,50,84],labels = ["cps", "sigma", "bg_cps"]):
    assert len(percentiles) == 3
    print("Best-fit parameters:")
    stats = {l:[] for l in labels}
    for i in range(len(labels)):
        mcmc = np.percentile(flat_samples[:, i], percentiles)
        q = np.diff(mcmc)
        print(f"{labels[i]}: {mcmc[1]:.3f} +{q[1]:.3f} -{q[0]:.3f}")
        stats[labels[i]] = (mcmc[1],q[0],q[1])
    return stats

def plot_mcmc_results(flat_samples,aperture_radii,flux,flux_err,
                      percentiles=[16,50,84],
                      labels = ["cps", "sigma", "bg_cps"],
                      r = np.arange(1.5,17.3,0.1)):

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

    plt.xlabel('aperture radius (as)')
    plt.ylabel('flux (cps)')
    plt.legend()
    plt.tight_layout()
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