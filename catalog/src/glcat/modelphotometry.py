import numpy as np
import pyarrow.parquet as parquet
import pyarrow as pa
from scipy.optimize import curve_fit
from scipy.stats import halfnorm
from collections import defaultdict


def model_photometry(
    aperture_radii: list[float],
    phot: dict[str, list[float]],
    which_band: str
) -> dict[str, list[float]]:
    """
    run model photometry, one band at a time
    """
    aperture_radii = np.asarray(aperture_radii)
    n_rows = len(next(iter(phot.values())))  # length of the first column of phot, assuming all same length
    flux_cols = [f'{which_band}_CPS_A{i}' for i in range(len(aperture_radii))]
    flux_err_cols = [f'{which_band}_CPS_ERR_A{i}' for i in range(len(aperture_radii))]
    for col in flux_cols + flux_err_cols:
        if col not in phot:
            raise ValueError(f"missing column: {col}, can't run model photometry")

    results = defaultdict(list) # empty

    for row in range(n_rows):
        flux = np.array([phot[col][row] for col in flux_cols], dtype=np.float64)
        flux_err = np.array([phot[col][row] for col in flux_err_cols], dtype=np.float64)

        # this shouldn't happen unless it's masked like a hotspot
        if (np.any(np.isnan(flux)) or np.any(np.isnan(flux_err)) or
            np.any(np.isinf(flux)) or np.any(np.isinf(flux_err))):
            # skip modeling on bad rows (this might break if they're pyarrow nulls?)
            # None later gets changed to pyarrow null when building the table?
            for j in range(len(aperture_radii)):
                results[f"{which_band}_MDL_SE_A{j}"].append(None)
                results[f"{which_band}_MDL_RSL_A{j}"].append(None)
            for name in ["CPS", "SIGMA", "BKG_CPS", "CPS_SE", "SIGMA_SE", "BKG_CPS_SE", "DET"]:
                results[f"{which_band}_MDL_{name}"].append(None)
            continue

        flux_err[flux_err == 0] = np.inf # have to do this to avoid divide by 0
        try:
            guess = [np.mean(flux),
                    5 / 2.355,
                    (flux[-1] - flux[-2]) / (np.pi * (aperture_radii[-1] ** 2 - aperture_radii[-2] ** 2))]
            bounds = ((0, 0, 0), (1000, 10, 100))
            cf_params, cf_cov = curve_fit(
                f=gaussian_flux_model_curvefit,
                xdata=aperture_radii,
                ydata=flux,
                p0=guess,
                sigma=flux_err,
                absolute_sigma=True,
                check_finite=False,
                method='trf',
                bounds=bounds)
            cf_param_errs = np.sqrt(np.diag(cf_cov))
            dists = [halfnorm.rvs(p, s, 3200) for p, s in zip(cf_params, cf_param_errs)]
            cf_sigmas = [np.std(gaussian_flux_model_curvefit(r, *dists)) for r in aperture_radii]
            modeled_cps = gaussian_flux_model_curvefit(aperture_radii, *cf_params)
            cf_residuals = modeled_cps - flux
            det = 1 - (np.var(cf_residuals) / np.var(flux))

            results[f"{which_band}_MDL_CPS"].append(cf_params[0])
            results[f"{which_band}_MDL_SIGMA"].append(cf_params[1])
            results[f"{which_band}_MDL_BKG_CPS"].append(cf_params[2])
            results[f"{which_band}_MDL_CPS_SE"].append(cf_param_errs[0])
            results[f"{which_band}_MDL_SIGMA_SE"].append(cf_param_errs[1])
            results[f"{which_band}_MDL_BKG_CPS_SE"].append(cf_param_errs[2])
            results[f"{which_band}_MDL_DET"].append(det)
            for j in range(len(aperture_radii)):
                results[f"{which_band}_MDL_SE_A{j}"].append(cf_sigmas[j])
                results[f"{which_band}_MDL_RSL_A{j}"].append(cf_residuals[j])

        except Exception as e:
            print(f"curve fit failed on row {row}: {e}")
            for j in range(len(aperture_radii)):
                results[f"{which_band}_MDL_SE_A{j}"].append(None)
                results[f"{which_band}_MDL_RSL_A{j}"].append(None)
            for name in ["CPS", "SIGMA", "BKG_CPS", "CPS_SE", "SIGMA_SE", "BKG_CPS_SE", "DET"]:
                results[f"{which_band}_MDL_{name}"].append(None)

    return results


def gaussian_flux_model_curvefit(r, total_flux, sigma, background):
    return (
        total_flux
        * (1 - np.exp(-r**2 / (2 * sigma**2)))
        + background * np.pi * r**2
    )
