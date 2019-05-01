import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.engine.url import URL
from sqlalchemy.pool import NullPool

_REGULAR_CHR = {"chr{}".format(x) for x in (list(range(1, 23)) + ["X", "Y"])}  # range(1,23) = 1,2,...,22


class GenomeBrowserClient:
    # If you are going to use the `local_hg19` configuration, make sure you have created such a user in MySQL:
        # For MySQL 5.6
        #   GRANT SELECT PRIVILEGES ON hg19.* To 'bud'@'localhost' IDENTIFIED BY 'earth';
        #   GRANT SELECT PRIVILEGES ON hgmd_pro.* To 'bud'@'localhost' IDENTIFIED BY 'earth';
        # For MySQL 5.7
        #   CREATE USER 'bud'@'localhost' IDENTIFIED BY 'earth';
        #   GRANT SELECT ON hg19.* TO 'bud'@'localhost';
        #   GRANT SELECT ON hgmd_pro.* TO 'bud'@'localhost';
        #   FLUSH PRIVILEGES;
    __db_url = dict(
        local_hg19=dict(
            drivername='mysql+pymysql',
            host='localhost',
            port='3306',
            username='bud',
            password='earth',
            database='hg19',
            query={'charset': 'utf8'}
        ),

        remote_hg19=dict(
            drivername='mysql+pymysql',
            host='genome-mysql.cse.ucsc.edu',
            port='3306',
            username='genome',
            password='',
            database='hg19',
            query={'charset': 'utf8'}
        ),
    )

    def __init__(self, config_key):
        # For `poolclass`, see http://stackoverflow.com/a/8705750
        self.db = create_engine(URL(**GenomeBrowserClient.__db_url[config_key]), poolclass=NullPool)
        self.conn = self.db.connect()

        # Subtraction between integer values, where one is of type UNSIGNED, produces an unsigned result by default.
        #   If the difference is negative, an error results because it must be unsigned.
        # Coordinates are unsigned int. We'll use subtraction between coordinates to get TSS distances,
        #   so we must enable this mode.
        self.conn.execute("SET sql_mode = 'NO_UNSIGNED_SUBTRACTION'")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()
        self.db.dispose()

    def select_tfbs(self, rsid):
        snps = ", ".join("'{}'".format(x) for x in rsid)
        chroms = ", ".join("'{}'".format(x) for x in _REGULAR_CHR)
        clazz = "'single'"

        query = ("SELECT s.name, GROUP_CONCAT(tf.name) as tfName "
                 "FROM snp146 as s "
                 "LEFT OUTER JOIN wgEncodeRegTfbsClusteredV3 as tf "
                 "ON tf.bin = s.bin "
                 "  AND s.chromStart BETWEEN tf.chromStart AND tf.chromEnd - 1 "
                 "  AND tf.chrom = s.chrom "
                 "WHERE s.name IN ({snps}) AND "
                 "  s.chrom IN ({chroms}) AND "
                 "  s.class = {clazz} "
                 "GROUP BY s.name".format(snps=snps, chroms=chroms, clazz=clazz))

        rows = self.conn.execute(query)

        df = pd.DataFrame(rows.fetchall())
        df.columns = rows.keys()

        return df


def binary_encode_tfbs(dfm, target_colname="tfName", value_sep=',', dest_colname_prefix=None):
    """
    Binary-encode categorical column `target_colname` separated by `value_sep`
    in data frame `dfm` to multiple binary columns with the same prefix `dest_colname_prefix`

    MySQL returns `GROUP_CONCAT(tf.name)` in `tfName` column in a comma-separated string.
        E.g. `ARID3A,ATF1,ATF2` stands for 3 TFs
    This function separates this string by commas and 3 new columns,
        `tf_ARID3A`, `tf_ATF1` and `tf_ATF2` would be 1 for this SNP

    :param dfm: the data frame
    :param target_colname: the name of the categorical column whose values would be encoded
    :param value_sep: the separator of the categorical values
    :param dest_colname_prefix: the prefix of the binary columns after one-hot encoding
    :return: the binary encoded dataframe
    """
    dummies = dfm.loc[:, target_colname].str.get_dummies(sep=value_sep)

    if dest_colname_prefix is not None:
        # Add a prefix to all column names
        dummies = dummies.add_prefix(dest_colname_prefix)

    dfm = pd.concat([dfm, dummies], axis=1).drop(target_colname, axis=1)

    return dfm
