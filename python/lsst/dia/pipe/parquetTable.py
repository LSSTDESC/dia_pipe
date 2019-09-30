import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


class ParquetTable(object):
    """Thin wrapper to pyarrow's ParquetFile object
    Call `toDataFrame` method to get a `pandas.DataFrame` object,
    Parameters
    ----------
    filename : str
        Path to Parquet file.
    """
    def __init__(self, filename=None, dataFrame=None):
        if filename is not None:
            self._pf = pq.ParquetFile(filename)
            self._df = None
        elif dataFrame is not None:
            self._df = dataFrame
            self._pf = None
        else:
            raise ValueError('Either filename or dataFrame must be passed.')

    def write(self, filename):
        """Write pandas dataframe to parquet
        Parameters
        ----------
        filename : str
            Path to which to write.
        """
        if self._df is None:
            raise ValueError('df property must be defined to write.')
        table = pa.Table.from_pandas(self._df)
        pq.write_table(table, filename, compression='none')

    def toDataFrame(self):
        """Get table as a pandas DataFrame
        """
        if self._df:
            return self._df

        return self._pf.read().to_pandas()
