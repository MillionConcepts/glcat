import os
import sys
import time
import numpy as np
import pyarrow.parquet as parquet
import pyarrow as pa
from scipy.optimize import curve_fit
from scipy.stats import halfnorm

def gaussian_flux_model_curvefit(r, total_flux, sigma, background):
    return (
        total_flux
        * (1 - np.exp(-r**2 / (2 * sigma**2)))
        + background * np.pi * r**2
    )

def model_fit_results(table, which_band, aperture_radii, bounds):
    """ iterate through rows (sources) in the catalog file and add model photometry stats """

    flux_cols = [f'{which_band}_CPS_A{i}' for i in range(len(aperture_radii))]
    flux_err_cols = [f'{which_band}_CPS_ERR_A{i}' for i in range(len(aperture_radii))]
    all_cols = flux_cols + flux_err_cols
    subset = table.select(all_cols)

    # New column names:
    # MDL_CPS, MDL_SIGMA, MDL_BKG_CPS = model fit parameters
    # MDL_CPS_SE, MDL_SIGMA_SE, MDL_BKG_CPS_SE = model parameters standard error
    # MDL_DET = coefficient of determination, quality metric
    # MDL_SE_Aix = standard error per aperture
    # MDL_RSL_Aix = residuals per aperture
    results = {
        f"{which_band}_MDL_CPS": [],
        f"{which_band}_MDL_SIGMA": [],
        f"{which_band}_MDL_BKG_CPS": [],
        f"{which_band}_MDL_CPS_SE": [],
        f"{which_band}_MDL_SIGMA_SE": [],
        f"{which_band}_MDL_BKG_CPS_SE": [],
        f"{which_band}_MDL_DET": [],
    }
    for i in range(len(aperture_radii)):
        results[f"{which_band}_MDL_SE_A{i}"] = []
        results[f"{which_band}_MDL_RSL_A{i}"] = []

    for i in range(subset.num_rows):
        row = subset.slice(i, 1).to_pandas().iloc[0]
        # we don't want to try and fit the model if flux or flux err are None or infinity
        # inf sometimes happens at hotspots
        if row.isnull().any() or np.isinf(row).any():
            for key in results:
                results[key].append(None)
            continue

        try:
            flux = np.array([row[c] for c in flux_cols])
            flux_err = np.array([row[c] for c in flux_err_cols])
            # avoid divide by zero in gaussian_flux_model_curvefit
            flux_err[flux_err == 0] = np.inf
            guess = [
                np.mean(flux),
                5 / 2.355,
                (flux[-1] - flux[-2]) / (np.pi * (aperture_radii[-1]**2 - aperture_radii[-2]**2))
            ]
            cf_params, cf_cov = curve_fit(
                f=gaussian_flux_model_curvefit,
                xdata=aperture_radii,
                ydata=flux,
                p0=guess,
                sigma=flux_err,
                absolute_sigma=True,
                check_finite=False,
                method='trf',
                bounds=bounds,
            )
            cf_param_errs = np.sqrt(np.diag(cf_cov))
            dists = [halfnorm.rvs(p, s, 3200) for (p, s) in zip(cf_params, cf_param_errs)]
            cf_sigmas = [
                np.std(gaussian_flux_model_curvefit(r, *dists))
                for r in aperture_radii
            ]
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

        except Exception:
            for key in results:
                results[key].append(None)
    return results

def run(inputs):
    eclipse, leg, band, input_dir, output_dir = inputs
    start = time.time()
    # update if aperture ranges change, could be input like in glcat?
    aperture_radii = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30])
    bounds = ((0, 0, 0), (1000, 10, 100))

    catalog_in = os.path.join(input_dir, f'e{eclipse}', f'e{eclipse}-{band}d-{leg}-catalog.parquet')
    if not os.path.exists(catalog_in):
        print(f"file missing for band={band}, leg={leg}, eclipse={eclipse}")
        print(catalog_in)
        return

    table = parquet.read_table(catalog_in)

    all_results = {}
    for which_band in ('NUV', 'FUV'):
        band_results = model_fit_results(table, which_band, aperture_radii, bounds)
        all_results.update(band_results)
    # even dicts with just keys eval as true
    if all_results:
        new_arrays = [pa.array(all_results[key], type=pa.float32()) for key in all_results]
        new_names = list(all_results.keys())
        table = pa.table(table.columns + new_arrays, names=table.schema.names + new_names)

    catalog_out = os.path.join(output_dir, f'e{eclipse}', f'e{eclipse}-{band}d-{leg}-catalog.parquet')
    os.makedirs(os.path.dirname(catalog_out), exist_ok=True)
    parquet.write_table(table, catalog_out)

    end = time.time()
    print(f"out:{catalog_out}, sources:{table.num_rows}, in {end - start:.2f} sec")

def main():
    if len(sys.argv) != 6:
        print("You need inputs silly!")
        print("modelphotom <eclipse:int> <leg:int> <band:str> <input_dir> <output_dir>")
        sys.exit(1) # don't give error traceback
    eclipse = int(sys.argv[1])
    # the catalog files are not padded with 0s in front of leg number
    leg = int(sys.argv[2])
    band = sys.argv[3].lower()
    input_dir = sys.argv[4]
    output_dir = sys.argv[5]
    run((eclipse, leg, band, input_dir, output_dir))
