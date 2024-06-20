import numpy as np
import emcee
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

ncat = pd.read_csv('ncat_240616.csv',index_col=None)

def gaussian_flux_model(theta, r):
    total_flux, sigma, background = theta
    #total_flux, background = theta
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

def fit_aperture_curve(aperture_radii,flux,flux_err,
                       showplots=False, # show the walker plots
                       showstats=False, # print summary results
                       nsteps = 5000, # number of MCMC steps
                       burnin = 200, # number of burn-in steps
                      ):
    # Set up the MCMC sampler
    ndim = 3  # Number of parameters to fit (total_flux, sigma, background)
    nwalkers = 32  # Number of MCMC walkers
    
    # Initialize the walkers
    initial_guess = [np.mean(flux), 5/2.355, (flux[-1]-flux[-2])/(np.pi * (aperture_radii[-1]**2-aperture_radii[-2]**2))]
    pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)
    
    # Run the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(aperture_radii, flux, flux_err))
    burn_in = sampler.run_mcmc(pos, burnin)
    sampler.reset()
    sampler.run_mcmc(burn_in, nsteps)

    if showplots:
        # Plot the MCMC results
        fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
        samples = sampler.get_chain()
        labels = ["Total Flux", "Sigma", "Background"]
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("Step number")
        plt.tight_layout()
        plt.show()

    # Print the best-fit parameters
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    if showstats:
        print("Best-fit parameters:")
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            print(f"{labels[i]}: {mcmc[1]:.3f} +{q[1]:.3f} -{q[0]:.3f}")
    return flat_samples

aperture_radii = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3])

np.random.seed(1965)
models = pd.DataFrame({})
for i in tqdm(np.arange(len(ncat))):
    flux = np.array([ncat[f'FUV_CPS_APER{a}'].iloc[i] for a,r in enumerate(aperture_radii)])
    flux_err = np.array([ncat[f'FUV_CPS_ERR_APER{a}'].iloc[i] for a,r in enumerate(aperture_radii)])
    ix = np.where(np.isfinite(flux))
    flat_samples = fit_aperture_curve(aperture_radii[ix],flux[ix],flux_err[ix],nsteps=200)

    a,b,c = (np.percentile(flat_samples[:, 0], [16, 50, 84])[1],
             np.percentile(flat_samples[:, 1], [16, 50, 84])[1],
             np.percentile(flat_samples[:, 2], [16, 50, 84])[1])
    
    q_a = np.diff(np.percentile(flat_samples[:, 0], [16, 50, 84]))
    q_b = np.diff(np.percentile(flat_samples[:, 1], [16, 50, 84]))
    q_c = np.diff(np.percentile(flat_samples[:, 2], [16, 50, 84]))
    
    models = pd.concat([models,
                        pd.DataFrame({'cps':[a],
                                      'cps_lower_conf':[q_a[0]],
                                      'cps_upper_conf':[q_a[1]],
                                      'sigma':[b],
                                      'sigma_lower_conf':[q_b[0]],
                                      'sigma_upper_conf':[q_b[1]],
                                      'bg_cps':[c],
                                      'bg_cps_lower_conf':[q_c[0]],
                                      'bg_cps_upper_conf':[q_c[1]],
                                      'bg_cps_annulus':[(flux[-1]-flux[-2])/(np.pi * (aperture_radii[-1]**2 - aperture_radii[-2]**2))],
                                      })],
                       axis=0)

#    plt.figure()
#    plt.errorbar(aperture_radii,flux,yerr=flux_err,fmt='k-',label='data')
#    plt.plot(aperture_radii,
#             gaussian_flux_model([a,b,c,], aperture_radii),'b:',
#             label=f'flux,sigma,bg= ({a:.3f} {b:.3f} {c:.3f})')
#    plt.legend()
#    plt.semilogy();

models.to_csv('ncat_fuv_fit_240616.csv',index=False)