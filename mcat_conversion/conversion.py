import json
from pathlib import Path

from hostess.aws.s3 import Bucket
import pandas as pd
import pdr
import pyarrow as pa
from pyarrow import parquet as pq

from column_settings import CASTCOLS, DICTCOLS


def mcat_to_pyarrow(mcat_path: Path) -> tuple[pa.Table, list[str]]:
    """
    Convert the MCAT .fits.gz file at `mcat_path` to a pyarrow table.
    Returns a tuple like (table, columns we tried to downcast but didn't).
    TODO: The max absolute error is currently hardcoded as a magic number
     (2.5e-7).
    """
    mcat = pdr.fastread(mcat_path)
    sources, newcols = mcat.GALEX_MERGED_SOURCE_LIST, {}
    warncols = []
    for colname in sources.columns:
        if colname in CASTCOLS:
            newcols[colname] = sources[colname].astype(CASTCOLS[colname][0])
            # TODO, maybe: add something to explicitly catch cast overflow
            #  errors...although the error check should always catch it anyway.
            if abs(sources[colname] - newcols[colname]).max() > 2.5e-7:
                warncols.append(colname)
                newcols[colname] = sources[colname]
        else:
            newcols[colname] = sources[colname].copy()
    del sources
    # TODO, maybe: reconstructing a dataframe is inefficient. I was just too
    #  lazy to type out all stuff to create pyarrow Field/Schema
    #  specifications.
    pa_tab = pa.Table.from_pandas(pd.DataFrame(newcols), preserve_index=False)
    jsons = {}
    for k, v in mcat.metadata.items():
        jsons[k] = json.dumps(dict(v)).encode('utf-8')
    return pa_tab.replace_schema_metadata(jsons), warncols


def convert_mcats(
    bucket_name: str,
    keys: list[str],
    temp_path: Path,
    outpre: str,
    logfile: Path
) -> dict[str, list[str]]:
    """
    Retrieve all MCAT .fits.gz files at keys `keys` in bucket `bucket_name`,
    write them to a temp path, convert them to parquet files, scratch them
    to disk, upload them to the same bucket under prefix `outpre`, then clean
    up all temp files. Performs simple success/failure logging to `logfile`.
    Returns a dict showing which columns in which files we didn't downcast as
    planned due to error (this is mostly for intermediate diagnosis and can
    probably be removed).

    Logging is super simple. CSV file, no headers, one line per event.
    Columns are: status, event name, exception text (empty on success).
    If the function fails before attempting to convert any individual file
    (e.g., while fetching from the bucket due to permissions errors), the
    event name will be 'initializing'.

    TODO: There is a deficiency in logging /
     handling: it only logs the failed file, then stops. It should continue
     through all other files instead, or at least log them all.
    """
    new_keys, logname, warncols = [], 'initializing', {}
    try:
        bucket = Bucket(bucket_name)
        bucket.get(keys, [temp_path / k for k in keys])
        for k in keys:
            name = k.replace("fits.gz", "parquet")
            pa_tab, wcols = mcat_to_pyarrow(temp_path / k)
            pq.write_table(
                pa_tab, temp_path / name, use_dictionary=DICTCOLS
            )
            warncols[k] = wcols
            new_keys.append(f"{outpre}/{name}")
        bucket.put(
            [temp_path / Path(k).name for k in new_keys], new_keys
        )
        with open(logfile, "a") as stream:
            for name in new_keys:
                stream.write(f"success,{Path(name).name},\n")
        return warncols
    except Exception as ex:
        with open(logfile, "a") as stream:
            stream.write(
                f"failure,{Path(name).name}"
                f",{str(ex).replace('\n','_').replace(',','_')}\n"
            )
        raise
    finally:
        for k in keys:
            if (tgz := (temp_path / k)).exists():
                tgz.unlink()
        for n in new_keys:
            if (tpq := (temp_path / Path(n).name)).exists():
                tpq.unlink()

