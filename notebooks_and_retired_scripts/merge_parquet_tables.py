"""Merge the records of UPDATE into TABLE.  The original version of TABLE
is renamed to {TABLE}.bak and the updated version replaces TABLE.
UPDATE must have exactly the same schema as TABLE."""

import argparse
import os
import sys
import tempfile

from pathlib import Path
from pyarrow import parquet, Table

SC = parquet.SortingColumn



def rewrite_table(path, table, *args, **kwargs):
    """Replace the file named "path" with a parquet file containing
       the contents of "table".  If "path" doesn't exist, create it,
       otherwise rename it with a ".bak" suffix; in that case, any
       existing file named "{path}.bak" will be overwritten."""

    target_dir = path.parent
    target_stem = path.stem

    backup_path = path.with_name(path.name + ".bak")

    # Infuriatingly, there is no way to create an Arrow NativeFile
    # directly from a file descriptor.  (You can create one from a
    # Python file object, but then you have multiple layers of buffering.)
    tempfd, temppath = tempfile.mkstemp(
        dir=target_dir,
        prefix=target_stem,
        suffix=".partial",
    )
    os.close(tempfd)
    temppath = Path(temppath)
    try:
        parquet.write_table(table, temppath, *args, **kwargs)
        try:
            path.rename(backup_path)
        except FileNotFoundError:
            pass

        temppath.rename(path)

    except:
        try:
            temppath.unlink()
        except FileNotFoundError:
            pass
        raise


def require_existing_file(f):
    try:
        r = f.resolve(strict=True)
    except FileNotFoundError:
        sys.stderr.write(f"{f}: error: No such file or directory\n")
        sys.exit(1)
    if not r.is_file():
        sys.stderr.write(f"{f}: error: Not a regular file\n")
        sys.exit(1)
    return r


def require_matching_schemas(schema1, schema2, fname1, fname2):
    if not schema1.equals(schema2, check_metadata=True):
        import difflib
        sys.stderr.write("error: schema mismatch:\n")
        s1 = schema1.to_string().splitlines()
        s2 = schema2.to_string().splitlines()
        n = max(len(s1), len(s2))
        sys.stderr.writelines(
            f"{line}\n"
            for line in difflib.unified_diff(
                s1, s2,
                str(fname1), str(fname2),
                n=n,
                lineterm="",
            )
            if not line.startswith("@@")
        )
        sys.exit(1)


def require_columns(needed, available, label):
    missing = [c for c in needed if c not in available]
    if missing:
        missing = ", ".join(missing)
        available = ", ".join(available)
        sys.stderr.write(
            f"error: {label} columns not found in tables: {missing}\n"
            f"note: available columns are: {available}\n"
        )
        sys.exit(1)


def parse_sort_order(spec):
    order = []
    spec = spec.strip()
    if spec == '':
        return order
    for colspec in spec.split(","):
        colspec = colspec.strip()
        if colspec.startswith("+"):
            order.append((colspec[1:], "ascending"))
        elif colspec.startswith("-"):
            order.append((colspec[1:], "descending"))
        else:
            raise ValueError(f"invalid --sort colspec {colspec!r}")
    return order


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table", type=Path, help="Parquet table to be updated.")
    ap.add_argument("update", type=Path, help="Table containing new records.")
    ap.add_argument("key_column", nargs="+",
                    help="One or more columns that, together, constitute"
                    " a primary key for the table being updated."
                    " All rows of the original, whose key columns match"
                    " at least one row of the update, will be replaced.")
    ap.add_argument("-z", "--compression",
                    help="Compression algorithm to use when writing the"
                    " updated table; passed directly to parquet.write_table."
                    " Valid values are 'none', 'snappy', 'gzip', 'brotli',"
                    " 'lz4', and 'zstd'; the default is 'snappy'.")
    ap.add_argument("--no-dictionary", dest="dictionary", action="store_false",
                    help="Disable dictionary encoding when writing"
                    " the updated table.")
    ap.add_argument("--sort", metavar="±col[,±col[,...]]", default="",
                    help="Sort the updated table by each of the specified"
                    " columns, in the sequence given.  Each must have a"
                    " leading + or -, specifying sorting in ascending or"
                    " descending order, respectively. If this option is"
                    " not given, the new table will not be sorted at all."
                    " Note that you must write --sort=-col,... if you want"
                    " to start with a column sorted in descending order."
                    " Due to a bug in the argparse library, --sort -col will"
                    " be treated as a syntax error.")

    args = ap.parse_args()
    try:
        sort_order = parse_sort_order(args.sort)
    except ValueError as e:
        ap.error(e.args[0])

    write_kwargs = {}
    if args.compression is not None:
        write_kwargs["compression"] = args.compression
    if not args.dictionary:
        write_kwargs["use_dictionary"] = False

    table_path = require_existing_file(args.table)
    update_path = require_existing_file(args.update)
    if table_path == update_path:
        sys.stderr.write(
            f"warning: {args.table} and {args.update} are the same file\n"
            f"note: it will be left unmodified\n"
        )
        sys.exit(0)

    with parquet.ParquetFile(table_path) as table, \
         parquet.ParquetFile(update_path) as update:
        table_schema = table.schema_arrow
        require_matching_schemas(
            table_schema,
            update.schema_arrow,
            args.table,
            args.update
        )

        key_columns = args.key_column
        table_columns = frozenset(table_schema.names)
        require_columns(key_columns, table_columns, "key")
        require_columns((o[0] for o in sort_order), table_columns, "sorting")

        # Read in the updates and construct the set of key tuples to be
        # filtered out of the original table.
        batches = []
        exclude_keys = set()
        for ub in update.iter_batches():
            exclude_keys.update(
                zip(*(ub.column(col) for col in key_columns))
            )
            batches.append(ub)

        # Read in the original table and filter out the rows to be
        # replaced.  Pyarrow expressions don't have any concept of
        # tuples as far as I can tell, so we have to construct a mask
        # array for the batch by hand.
        for tb in table.iter_batches():
            mask = [
                (key not in exclude_keys)
                for key in zip(*(tb.column(col) for col in key_columns))
            ]
            batches.append(tb.filter(mask))

    new_table = Table.from_batches(batches)

    if sort_order:
        new_table = new_table.sort_by(sort_order)
        write_kwargs["sorting_columns"] = SC.from_ordering(
            new_table.schema,
            sort_order,
        )

    rewrite_table(table_path, new_table, **write_kwargs)


main()
