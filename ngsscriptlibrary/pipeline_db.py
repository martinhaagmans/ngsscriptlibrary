"""Database IO functions for pipeline."""
import os
import time
import json
import sqlite3
from collections import namedtuple
from collections import defaultdict

import pandas as pd


class SampleSheet:
    """Create a samplesheet to run the MiSeq and the analysis pipeline."""

    def __init__(self, nullijst, serie, samplesheet, faciliteit=False):
        """Parse serie specific info."""
        self.nullijst = nullijst
        self.serie = serie
        self.samplesheet = samplesheet
        self.basedir = os.path.dirname(os.path.abspath(__file__))
        if faciliteit:
            self.barcodedict = {k: k for k in range(1, 193)}
            self.faciliteit = True
        if not faciliteit:
            self.barcodedict = self.create_barcode_hussle_dict()
            self.faciliteit = False
        self.adapterdict = self.create_adapter_dict()
        self.adaptercodes = self.create_adapter_code_dict()

    def create_barcode_hussle_dict(self):
        """Turn TSV-file to dict.

        Read tab seperated file with our barcodenumber and the official one.
        Return dict with our number as key and the official as value.
        """
        barcodedict = dict()
        with open(os.path.join(self.basedir, 'docs/barcodehussle.txt')) as f:
            for line in f:
                if not line.split():
                    continue
                else:
                    us, them = line.strip().split()
                    barcodedict[int(us)] = int(them)
        return barcodedict

    def create_adapter_dict(self):
        """Turn TSV-file to dict.

        Read tab seperated file with barcode, seqI7 and seqI5.
        Return dict with barcode as key and a tuple of sequences as value.
        """
        adapterdict = dict()

        with open(os.path.join(self.basedir, 'docs/adapters.txt')) as f:
            for line in f:
                if not line.split():
                    continue
                else:
                    index, left, right = line.strip().split('\t')
                    adapterdict[int(index)] = [left, right]
        return adapterdict

    def create_adapter_code_dict(self):
        """Turn TSV-file to dict.

        Read tab seperated file with numbered seqI7 and seqI5.
        Return dict with sequence as key and number as value.
        """
        adapterdict = dict()
        with open(os.path.join(self.basedir, 'docs/adaptercodes.txt')) as f:
            for line in f:
                if not line.split():
                    continue
                else:
                    code, sequence = line.strip().split('\t')
                    adapterdict[sequence] = code
        return adapterdict

    def parse_nullijst(self):
        """Parse file with sample, barcode, analysis return todo dictionary.

        Parse file and convert barcode from our number to seqi5 and seqi7 with
        Illumina numbers. Return dict with sampleID as key and all
        barcodenumbers/sequences and the diagnostic test requested as values.
        """
        out = dict()
        samples_parsed = list()

        with open(self.nullijst, 'r') as f:
            for line in f:
                sample, genesis, barcode = line.split()
                if sample not in samples_parsed:
                    samples_parsed.append(sample)
                elif sample in samples_parsed:
                    samples_parsed.append(sample)
                    sample = sample + 'A' * (samples_parsed.count(sample) - 1)
                barcode = int(barcode)
                barcode = self.barcodedict[barcode]
                genesis = genesis.replace('NGS-', '')
                out[sample] = dict()
                out[sample]['test'] = genesis
                out[sample]['barcode'] = dict()
                out[sample]['barcode']['i7code'] = self.adaptercodes[self.adapterdict[barcode][0]]
                out[sample]['barcode']['i7seq'] = self.adapterdict[barcode][0]
                out[sample]['barcode']['i5code'] = self.adaptercodes[self.adapterdict[barcode][1]]
                out[sample]['barcode']['i5seq'] = self.adapterdict[barcode][1]
        return out

    def writelist(self, list_to_print, f):
        """Write items in list to fileobject. If item > 1: print csv tuple."""
        for val in list_to_print:
            if not isinstance(val, int):
                if '[' in val and ']' in val:
                    f.write(str(val))
                elif len(val) == 2:
                    if not val[1] is None:
                        f.write('{},'.format(val[0]))
                        f.write(str(val[1]))
            elif isinstance(val, int):
                f.write('{}'.format(val))
            f.write('\n')

    def get_adapter_seq(self, adapterset):
        """Return adaptersequences for requested adapterset."""
        if adapterset == 'BIOO set 2':
            adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
            adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
            return adapter1, adapter2
        else:
            return ['Unknown', 'Unknown']

    def write_files(self, adapterset=None, readlength=None, analist=None):
        """Create samplesheet."""
        if analist is None:
            analist = 'Robert'
        if adapterset is None:
            adapterset = 'BIOO set 2'
        if readlength is None:
            readlength = 151
        date = str(time.strftime("%m-%d-%Y"))
        todo_nullijst = self.parse_nullijst()

        if self.faciliteit is False:
            captures = '_'.join(set([todo_nullijst[s]['test'] for s in todo_nullijst.keys()]))
        elif self.faciliteit is True:
            captures = 'Nvt'
            analist = 'LVGA'

        # Define samplesheet
        header = ['[Header]',
                  ['IEMFileVersion', 4],
                  ['Investigator Name', analist],
                  ['Experiment Name', self.serie],
                  ['Date', date],
                  ['Workflow', 'GenerateFASTQ'],
                  ['Application', 'FASTQ Only'],
                  ['Assay', adapterset],
                  ['Description', captures],
                  ['Chemistry', 'Amplicon']]

        reads = ['[Reads]', readlength, readlength]

        settings = ['[Settings]',
                    ['ReverseComplement', 0],
                    ['Adapter', self.get_adapter_seq(adapterset)[0]],
                    ['AdapterRead2', self.get_adapter_seq(adapterset)[1]]]

        header_samples = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                          'I7_Index_ID', 'index', 'I5_Index_ID', 'index2',
                          'Sample_Project', 'Description']

        with open(self.samplesheet, 'w') as f:
            self.writelist(header, f)
            f.write('\n')
            self.writelist(reads, f)
            f.write('\n')
            self.writelist(settings, f)
            f.write('\n')
            f.write('[Data]\n')
            f.write(','.join(header_samples))
            f.write('\n')
            for sample in todo_nullijst.keys():
                sample_out = ('{s},{s},,,{i7},{i7seq},{i5},{i5seq},{p},{t}\n'
                              .format(s=sample,
                                      i7=todo_nullijst[sample]['barcode']['i7code'],
                                      i7seq=todo_nullijst[sample]['barcode']['i7seq'],
                                      i5=todo_nullijst[sample]['barcode']['i5code'],
                                      i5seq=todo_nullijst[sample]['barcode']['i5seq'],
                                      p=self.serie, t=todo_nullijst[sample]['test']))
                f.write(sample_out)


class TargetDBCreator:
    """Create database tables for pipeline functionality.

    Every diagnostic test has a unique code that is linked to analyses to
    perform and target files to be used in given analyses.
    """

    def __init__(self, db):
        """Establish a connection with a database."""
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()

    def _execute_and_commit(self, statement):
        try:
            self.c.execute(statement)
        except sqlite3.OperationalError as e:
            print(e)
        else:
            self.conn.commit()

    def create_genesis_table(self):
        """Create table that links genesis code to analyses and files."""
        sql = """CREATE TABLE genesis (
        genesis TEXT NOT NULL,
        capture TEXT NOT NULL,
        pakket TEXT NOT NULL,
        panel TEXT,
        cnvscreening BOOLEAN NOT NULL DEFAULT 0,
        cnvdiagnostiek BOOLEAN NOT NULL DEFAULT 0,
        mozaiekdiagnostiek BOOLEAN NOT NULL DEFAULT 0,
        PRIMARY KEY (genesis))
        """
        self._execute_and_commit(sql)

    def create_aandoeningen_table(self):
        """Create table that links genesis code to phenotype."""
        sql = """CREATE TABLE aandoeningen (
        genesis TEXT NOT NULL,
        aandoening TEXT NOT NULL,
        PRIMARY KEY (genesis))
        """
        self._execute_and_commit(sql)

    def create_captures_table(self):
        """Create table for genes and size per capture."""
        sql = """CREATE TABLE captures (
        capture TEXT NOT NULL,
        versie INTEGER NOT NULL,
        oid INTEGER NOT NULL,
        verdund BOOLEAN NOT NULL,
        grootte INTEGER NOT NULL,
        genen TEXT,
        PRIMARY KEY (capture, versie))
        """
        self._execute_and_commit(sql)

    def create_pakketten_table(self):
        """Create table for genes and size per pakket."""
        sql = """CREATE TABLE pakketten (
        pakket TEXT NOT NULL,
        versie INTEGER NOT NULL,
        grootte INTEGER,
        genen TEXT,
        PRIMARY KEY (pakket, versie))
        """
        self._execute_and_commit(sql)

    def create_panels_table(self):
        """Create table for genes and size per panel."""
        sql = """CREATE TABLE panels (
        panel TEXT NOT NULL,
        versie INTEGER NOT NULL,
        grootte INTEGER,
        genen TEXT,
        PRIMARY KEY (panel, versie))
        """
        self._execute_and_commit(sql)

    def create_all(self):
        """Conveniently create all tables."""
        to_create = [self.create_genesis_table(),
                     self.create_aandoeningen_table(),
                     self.create_captures_table(),
                     self.create_pakketten_table(),
                     self.create_panels_table()]
        for func in to_create:
                func


class MetricsDBcreator:
    """Create database tables for pipeline output.

    Every diagnostic test has a unique code that is linked to analyses to
    perform and target files to be used in given analyses.
    """

    def __init__(self, db):
        """Establish a connection with a database."""
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()

    def _execute_and_commit(self, statement):
        try:
            self.c.execute(statement)
        except sqlite3.OperationalError as e:
            print(e)
        else:
            self.conn.commit()

    def create_samplesheet_table(self):
        """Create table to store sample sheet."""
        sql = """CREATE TABLE todo (
        SERIE TEXT NOT NULL,
        SAMPLE TEXT NOT NULL,
        genesis TEXT NOT NULL,
        capture TEXT NOT NULL,
        pakket TEXT NOT NULL,
        panel TEXT NOT NULL,
        cnvscreening INTEGER NOT NULL,
        cnvdiagnostiek INTEGER NOT NULL,
        mozaiekdiagnostiek INTEGER NOT NULL,
        PRIMARY KEY(SERIE, SAMPLE))
        """
        self._execute_and_commit(sql)

    def create_sanger_table(self):
        """Create table to store sanger fragments."""
        sql = """CREATE TABLE sangers (
        SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        TARGET TEXT NOT NULL,
        DATA TEXT NOT NULL,
        PRIMARY KEY(SAMPLE, SERIE, TARGET))
        """
        self._execute_and_commit(sql)

    def create_snpcheck_table(self):
        """Create table to store snpcheck."""
        sql = """CREATE TABLE snpcheck (
        SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        DATA TEXT NOT NULL,
        PRIMARY KEY(SAMPLE, SERIE))
        """
        self._execute_and_commit(sql)

    def create_alignmentmetrics_table(self):
        """Create table to store alignent metrics."""
        sql = """CREATE TABLE alignment (
        SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        CATEGORY TEXT NOT NULL,
        TOTAL_READS INTEGER NOT NULL,
        PF_READS INTEGER NOT NULL,
        PCT_PF_READS INTEGER NOT NULL,
        PF_NOISE_READS INTEGER NOT NULL,
        PF_READS_ALIGNED INTEGER NOT NULL,
        PCT_PF_READS_ALIGNED REAL NOT NULL,
        PF_ALIGNED_BASES INTEGER NOT NULL,
        PF_HQ_ALIGNED_READS INTEGER NOT NULL,
        PF_HQ_ALIGNED_BASES INTEGER NOT NULL,
        PF_HQ_ALIGNED_Q20_BASES INTEGER NOT NULL,
        PF_HQ_MEDIAN_MISMATCHES INTEGER NOT NULL,
        PF_MISMATCH_RATE REAL NOT NULL,
        PF_HQ_ERROR_RATE REAL NOT NULL,
        PF_INDEL_RATE REAL NOT NULL,
        MEAN_READ_LENGTH INTEGER NOT NULL,
        READS_ALIGNED_IN_PAIRS INTEGER NOT NULL,
        PCT_READS_ALIGNED_IN_PAIRS REAL NOT NULL,
        BAD_CYCLES INTEGER NOT NULL,
        STRAND_BALANCE REAL NOT NULL,
        PCT_CHIMERAS REAL NOT NULL,
        PCT_ADAPTER INTEGER NOT NULL,
        LIBRARY TEXT,
        READ_GROUP TEXT,
        PRIMARY KEY(SAMPLE,SERIE,CATEGORY))
        """
        self._execute_and_commit(sql)

    def create_hsmetrics_table(self):
        """Create table to store Hybrid Selection (HS) metrics."""
        sql = """CREATE TABLE hsmetrics (
        SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        BAIT_SET TEXT NOT NULL,
        GENOME_SIZE INTEGER NOT NULL,
        BAIT_TERRITORY INTEGER NOT NULL,
        TARGET_TERRITORY INTEGER NOT NULL,
        BAIT_DESIGN_EFFICIENCY REAL NOT NULL,
        TOTAL_READS INTEGER NOT NULL,
        PF_READS INTEGER NOT NULL,
        PF_UNIQUE_READS INTEGER NOT NULL,
        PCT_PF_READS INTEGER NOT NULL,
        PCT_PF_UQ_READS INTEGER NOT NULL,
        PF_UQ_READS_ALIGNED INTEGER NOT NULL,
        PCT_PF_UQ_READS_ALIGNED REAL NOT NULL,
        PF_BASES_ALIGNED INTEGER NOT NULL,
        PF_UQ_BASES_ALIGNED INTEGER NOT NULL,
        ON_BAIT_BASES INTEGER NOT NULL,
        NEAR_BAIT_BASES INTEGER NOT NULL,
        OFF_BAIT_BASES INTEGER NOT NULL,
        ON_TARGET_BASES INTEGER NOT NULL,
        PCT_SELECTED_BASES REAL NOT NULL,
        PCT_OFF_BAIT REAL NOT NULL,
        ON_BAIT_VS_SELECTED REAL NOT NULL,
        MEAN_BAIT_COVERAGE REAL NOT NULL,
        MEAN_TARGET_COVERAGE REAL NOT NULL,
        MEDIAN_TARGET_COVERAGE INTEGER NOT NULL,
        PCT_USABLE_BASES_ON_BAIT REAL NOT NULL,
        PCT_USABLE_BASES_ON_TARGET REAL NOT NULL,
        FOLD_ENRICHMENT REAL NOT NULL,
        ZERO_CVG_TARGETS_PCT REAL NOT NULL,
        PCT_EXC_DUPE INTEGER NOT NULL,
        PCT_EXC_MAPQ REAL NOT NULL,
        PCT_EXC_BASEQ INTEGER NOT NULL,
        PCT_EXC_OVERLAP INTEGER NOT NULL,
        PCT_EXC_OFF_TARGET REAL NOT NULL,
        FOLD_80_BASE_PENALTY TEXT NOT NULL,
        PCT_TARGET_BASES_1X REAL NOT NULL,
        PCT_TARGET_BASES_2X INTEGER NOT NULL,
        PCT_TARGET_BASES_10X INTEGER NOT NULL,
        PCT_TARGET_BASES_20X INTEGER NOT NULL,
        PCT_TARGET_BASES_30X INTEGER NOT NULL,
        PCT_TARGET_BASES_40X INTEGER NOT NULL,
        PCT_TARGET_BASES_50X INTEGER NOT NULL,
        PCT_TARGET_BASES_100X INTEGER NOT NULL,
        HS_LIBRARY_SIZE REAL,
        HS_PENALTY_10X INTEGER NOT NULL,
        HS_PENALTY_20X INTEGER NOT NULL,
        HS_PENALTY_30X INTEGER NOT NULL,
        HS_PENALTY_40X INTEGER NOT NULL,
        HS_PENALTY_50X INTEGER NOT NULL,
        HS_PENALTY_100X INTEGER NOT NULL,
        AT_DROPOUT INTEGER NOT NULL,
        GC_DROPOUT REAL NOT NULL,
        HET_SNP_SENSITIVITY REAL NOT NULL,
        HET_SNP_Q INTEGER NOT NULL,
        LIBRARY TEXT,
        READ_GROUP TEXT,
        PRIMARY KEY(SAMPLE,SERIE))
        """
        self._execute_and_commit(sql)

    def create_insertsizemetrics_table(self):
        """Create table to store insert size metrics."""
        sql = """CREATE TABLE insertsize (
        SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        MEDIAN_INSERT_SIZE INTEGER NOT NULL,
        MEDIAN_ABSOLUTE_DEVIATION INTEGER NOT NULL,
        MIN_INSERT_SIZE INTEGER NOT NULL,
        MAX_INSERT_SIZE INTEGER NOT NULL,
        MEAN_INSERT_SIZE REAL NOT NULL,
        STANDARD_DEVIATION REAL NOT NULL,
        READ_PAIRS INTEGER NOT NULL,
        PAIR_ORIENTATION TEXT NOT NULL,
        WIDTH_OF_10_PERCENT INTEGER NOT NULL,
        WIDTH_OF_20_PERCENT INTEGER NOT NULL,
        WIDTH_OF_30_PERCENT INTEGER NOT NULL,
        WIDTH_OF_40_PERCENT INTEGER NOT NULL,
        WIDTH_OF_50_PERCENT INTEGER NOT NULL,
        WIDTH_OF_60_PERCENT INTEGER NOT NULL,
        WIDTH_OF_70_PERCENT INTEGER NOT NULL,
        WIDTH_OF_80_PERCENT INTEGER NOT NULL,
        WIDTH_OF_90_PERCENT INTEGER NOT NULL,
        WIDTH_OF_99_PERCENT INTEGER NOT NULL,
        LIBRARY TEXT,
        READ_GROUP TEXT,
        PRIMARY KEY(SAMPLE,SERIE))
        """
        self._execute_and_commit(sql)

    def create_callable_summary_table(self):
        """Create table to store Callable Loci report."""
        sql = """CREATE TABLE callable
        (REF_N INTEGER NOT NULL,
        CALLABLE INTEGER NOT NULL,
        NO_COVERAGE INTEGER NOT NULL,
        LOW_COVERAGE INTEGER NOT NULL,
        EXCESSIVE_COVERAGE INTEGER NOT NULL,
        POOR_MAPPING_QUALITY INTEGER NOT NULL,
        TARGET TEXT NOT NULL,
        SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        PRIMARY KEY(SAMPLE, TARGET, SERIE))
        """
        self._execute_and_commit(sql)

    def create_perc_covered_table(self):
        """Create table to store percentage of target covered >30x."""
        sql = """CREATE TABLE procenttargetcovered
        (SAMPLE TEXT NOT NULL,
        TARGET TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        PERCENTAGE REAL NOT NULL,
        PRIMARY KEY(SAMPLE, TARGET, SERIE))
        """
        self._execute_and_commit(sql)

    def create_mean_std_coverage_table(self):
        """Create table to store informative mean and std dev."""
        sql = """CREATE TABLE mean
        (SAMPLE TEXT NOT NULL,
        TARGET TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        MEAN INTEGER NOT NULL,
        STDEV INTEGER NOT NULL,
        PRIMARY KEY(SAMPLE, TARGET, SERIE))
        """
        self._execute_and_commit(sql)        

    def create_base_perc_reads(self):
        """Create table to store per base percentages of R1."""
        sql = """CREATE TABLE basepercentages
        (SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        A REAL NOT NULL,
        T REAL NOT NULL,
        C REAL NOT NULL,
        G REAL NOT NULL,
        PRIMARY KEY(SAMPLE, SERIE))
        """
        self._execute_and_commit(sql)

    def create_risk_score_table(self):
        """Create table to store risk scores."""
        sql = """CREATE TABLE riskscore
        (SAMPLE TEXT NOT NULL,
        SERIE TEXT NOT NULL,
        TARGET TEXT NOT NULL,
        score real NOT NULL,
        genotypes TEXT NOT NULL,
        PRIMARY KEY(SAMPLE, SERIE))
        """
        self._execute_and_commit(sql)        

    def create_all(self):
        to_create = [self.create_samplesheet_table(),
                     self.create_sanger_table(),
                     self.create_alignmentmetrics_table(),
                     self.create_hsmetrics_table(),
                     self.create_insertsizemetrics_table(),
                     self.create_callable_summary_table(),
                     self.create_perc_covered_table(),
                     self.create_base_perc_reads(),
                     self.create_snpcheck_table(),
                     self.create_mean_std_coverage_table(),
                     self.create_risk_score_table()]

        for func in to_create:
                func


class MetricsDBReader:

    def __init__(self, db, sampleID, serie, target):
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()
        self.sampleID = sampleID
        self.serie = serie
        self.target = target

    def get_alignmetrics(self):
        read1 = '''SELECT DISTINCT
        TOTAL_READS,
        round(PCT_CHIMERAS*100, 1),
        round(MEAN_READ_LENGTH, 1),
        round(PCT_READS_ALIGNED_IN_PAIRS*100, 1)
        FROM {table}
        WHERE (CATEGORY = 'FIRST_OF_PAIR'
        AND sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='alignment', dnr=self.sampleID, serie=self.serie)

        read2 = '''SELECT DISTINCT
        round(MEAN_READ_LENGTH, 1),
        round(PCT_READS_ALIGNED_IN_PAIRS*100, 1)
        FROM {table}
        WHERE (CATEGORY = 'SECOND_OF_PAIR'
        AND sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='alignment', dnr=self.sampleID, serie=self.serie)

        self.c.execute(read1)
        read1out = [val for tup in self.c.fetchall() for val in tup]

        self.c.execute(read2)
        read2out = [val for tup in self.c.fetchall() for val in tup]
        return read1out + read2out

    def get_hsmetrics(self):
        sql = '''SELECT DISTINCT
        TARGET_TERRITORY,
        round(PCT_PF_UQ_READS*100, 1),
        round(PCT_SELECTED_BASES*100, 1)
        FROM {table}
        WHERE (sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='hsmetrics', dnr=self.sampleID, serie=self.serie)
        self.c.execute(sql)
        out = [val for tup in self.c.fetchall() for val in tup]
        return out

    def get_alignmetrics_for_seriereport(self):
        read1 = '''SELECT DISTINCT
        round(MEAN_READ_LENGTH)
        FROM {table}
        WHERE (CATEGORY = 'FIRST_OF_PAIR'
        AND sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='alignment', dnr=self.sampleID, serie=self.serie)

        read2 = '''SELECT DISTINCT
        round(MEAN_READ_LENGTH)
        FROM {table}
        WHERE (CATEGORY = 'SECOND_OF_PAIR'
        AND sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='alignment', dnr=self.sampleID, serie=self.serie)

        self.c.execute(read1)
        read1out = [val for tup in self.c.fetchall() for val in tup]

        self.c.execute(read2)
        read2out = [val for tup in self.c.fetchall() for val in tup]
        return read1out + read2out

    def get_hsmetrics_for_seriereport(self):
        sql = '''SELECT DISTINCT
        round(PCT_PF_UQ_READS*100),
        round(PCT_SELECTED_BASES*100),
        round(MEAN_TARGET_COVERAGE)
        FROM {table}
        WHERE (sample = '{dnr}'
        AND serie = '{serie}')
        '''.format(table='hsmetrics', dnr=self.sampleID, serie=self.serie)
        self.c.execute(sql)
        out = [val for tup in self.c.fetchall() for val in tup]
        return out

    def get_callable(self):
        total_sql = """SELECT REF_N, CALLABLE,
        NO_COVERAGE, LOW_COVERAGE, EXCESSIVE_COVERAGE
        FROM callable
        WHERE (SAMPLE='{s}' AND TARGET='{t}')
        """.format(s=self.sampleID, t=self.target)

        callable_sql = """SELECT CALLABLE
        FROM callable
        WHERE (SAMPLE='{s}' AND TARGET='{t}')
        """.format(s=self.sampleID, t=self.target)

        self.c.execute(total_sql)
        all_out = [val for tup in self.c.fetchall() for val in tup]

        self.c.execute(callable_sql)
        callable_out = [val for tup in self.c.fetchall() for val in tup]

        return (sum(callable_out) / sum(all_out)) * 100

    def get_perc_target_covered(self, diagnostic_target=None):
        if diagnostic_target is None:
            diagnostic_target = self.target
        sql = """SELECT PERCENTAGE
        FROM procenttargetcovered
        WHERE (SAMPLE='{sa}' AND SERIE='{se}' AND TARGET='{t}')
        """.format(sa=self.sampleID, se=self.serie, t=diagnostic_target)
        self.c.execute(sql)
        return self.c.fetchone()[0]

    def get_sanger_fragments(self, diagnostic_target=None):
        if diagnostic_target is None:
            diagnostic_target = self.target
        sql = """SELECT DATA
        FROM sangers
        WHERE (SAMPLE='{sa}' AND SERIE='{se}' AND TARGET='{t}')
        """.format(sa=self.sampleID, se=self.serie, t=diagnostic_target)
        self.c.execute(sql)
        sangers = json.loads(self.c.fetchone()[0])
        return sangers

    def get_snpcheck(self):
        sql = """SELECT DATA
        FROM snpcheck
        WHERE (SERIE='{}' AND SAMPLE='{}')
        """.format(self.serie, self.sampleID)
        self.c.execute(sql)
        snpcheck = json.loads(self.c.fetchone()[0])
        return snpcheck

    def get_mean_and_stdev(self):
        sql = """SELECT MEAN, STDEV
        FROM mean
        WHERE (SAMPLE='{sa}' AND SERIE='{se}' AND TARGET='{t}')
        """.format(sa=self.sampleID, se=self.serie, t=self.target)
        self.c.execute(sql)
        return self.c.fetchone()
    
    def get_risk_score(self):
        sql = """SELECT score
        FROM riskscore
        WHERE (SAMPLE='{sa}' AND SERIE='{se}' AND TARGET='{t}')
        """.format(sa=self.sampleID, se=self.serie, t=self.target)
        self.c.execute(sql)
        return self.c.fetchone()[0]        

    def get_all(self):
        print('Alignment: {}'.format(self.get_alignmetrics()),
              'HSmetrics: {}'.format(self.get_hsmetrics()),
              'Callable: {}'.format(self.get_callable()),
              'Perctargetcovered: {}'.format(self.get_perc_target_covered())
              )


def perc_target_covered2db(sampleID, percentage, pakket, serie, db):
    "Create database connection and write perc. target covered >30x to table."
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = '''INSERT INTO procenttargetcovered
    (SAMPLE, TARGET, SERIE, PERCENTAGE)
    VALUES ('{}', '{}', '{}', '{}')
    '''.format(sampleID, pakket, serie, percentage)
    try:
        c.execute(sql)
    except sqlite3.IntegrityError:
        pass
    else:
        conn.commit()
    conn.close()


def summary2db(sampleID, data, target, serie, db):
    "Insert dict with CallableLoci summary into database."
    data['TARGET'] = target
    data['SAMPLE'] = sampleID
    data['SERIE'] = serie
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = """INSERT INTO callable ({headers})
    VALUES ('{values}');
    """.format(headers=', '.join(data.keys()),
               values="""', '""".join(data.values()))
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        print(e)
    else:
        conn.commit()
    conn.close()


def metrics2db(df, db, table):
    "Create database connection and write picard results to table."
    conn = sqlite3.connect(db)
    try:
        df.to_sql(table, conn, if_exists='append')
    except sqlite3.OperationalError as e:
        print(e)
    conn.close()


def samplesheetinfo2db(info, sample, db):
    sql = '''INSERT INTO todo
    VALUES ('{}', '{}', '{}', '{}','{}', '{}', {}, {}, {})
    '''.format(info['serie'], sample, info['genesis'], info['capture'],
               info['pakket'], info['panel'], int(info['cnvscreening']),
               int(info['cnvdiagnostiek']), int(info['mozaiek']))
    conn = sqlite3.connect(db)
    c = conn.cursor()
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        print(e)
    else:
        conn.commit()
    conn.close()


def sangers2db(data, sample, serie, target, db):
    sql = '''INSERT INTO sangers
    VALUES ('{}', '{}', '{}', '{}')
    '''.format(serie, sample, target, json.dumps(data))
    conn = sqlite3.connect(db)
    c = conn.cursor()
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        print(e)
    else:
        conn.commit()
    conn.close()


def sample_in_db(sample, serie, db, table):
    """Check if sample-serie combi is in database. Return boolean."""
    sql = """SELECT * FROM {}
    WHERE(SAMPLE='{}' AND SERIE='{}')
    """.format(table, sample, str(serie))

    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(sql)

    if c.fetchone():
        conn.close()
        return True
    else:
        conn.close()
        return False


def base_perc_reads2db(sampleID, serie, percentages, db):
    "Create database connection and write basepercentages to table."
    a, t, c, g = percentages
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = '''INSERT INTO basepercentages
    (SAMPLE, SERIE, A, T, C, G)
    VALUES ('{}', '{}', {}, {}, {}, {})
    '''.format(sampleID, serie, a, t, c, g)
    try:
        cursor.execute(sql)
    except sqlite3.IntegrityError:
        pass
    else:
        conn.commit()
    conn.close()


def get_baseperc_reads_serie(serie, db):
    "Get basepercentages reads for serie. Return list of namedtuples."
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = '''SELECT SAMPLE, A, T, C, G 
    FROM basepercentages
    WHERE SERIE='{}'
    '''.format(serie)
    c.execute(sql)
    db_out = c.fetchall()
    conn.close()
    out = list()
    BaseP = namedtuple("BaseP", ['sample', 'a', 't', 'c', 'g'])
    for _ in db_out:
        sample, a, t, c, g = _
        __ = BaseP(sample, a, t, c, g)
        out.append(__)

    return out


def get_df_baseperc_reads_serie(serie, db):
    "Get basepercentages reads for serie. Return pandas df."
    sql = '''SELECT SAMPLE, A, T, C, G 
    FROM basepercentages
    WHERE SERIE='{}'
    '''.format(serie)
    df = pd.read_sql(sql, con=sqlite3.connect(db), index_col='SAMPLE')
    return df


def get_sangers_serie_from_db(serie, db=None):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = 'SELECT SAMPLE, TARGET, DATA FROM sangers WHERE SERIE="{}"'.format(str(serie))
    c.execute(sql)
    return c.fetchall()


def sangerdb_to_dataframe(data, sample):
    sanger_list = list()

    for __ in data:
        gene, chromosome, sanger_start, sanger_end, target_start, target_end = __
        sanger_list.append({'gen': gene, 
                            'chr': chromosome, 
                            'sanger_start': sanger_start, 
                            'sanger_end': sanger_end, 
                            'target_start': target_start, 
                            'target_end': target_end})

    df = pd.DataFrame(sanger_list)
    df['SAMPLE'] = sample
    return df


def combine_all_sangerdb_to_df(dbout, max_sangers=None):

    if max_sangers is None:
        max_sangers = 50

    dflist = list()
    failed = list()

    for _ in dbout:
        sample, target, data, *_ = _
        data = json.loads(data)

        if 'Geen sangers:' in data:
            continue
        elif len(data) >= max_sangers:
            failed.append(sample)
        else:
            dflist.append(sangerdb_to_dataframe(data, '{}:{}'.format(sample, target)))
 
    if dflist:
        df = pd.concat(dflist)
        df.set_index(['chr', 'target_start', 'target_end'], inplace=True)
        df.sort_index(inplace=True)

    elif not dflist:
        df = pd.DataFrame()
    
    return df, failed


def get_sangerfrag_min_max(df):
    df['max'], df['min'] = float('NaN'), float('NaN')
    for target in df.index:
        try:
            df.loc[target, 'max'] = df.loc[target]['sanger_end'].transpose().max()
        except AttributeError:
            df.loc[target, 'max'] = df.loc[target]['sanger_end']
        try:
            df.loc[target, 'min'] = df.loc[target]['sanger_start'].transpose().min()
        except AttributeError:            
            df.loc[target, 'min'] = df.loc[target]['sanger_start']            
    return df


def group_samples_by_sanger(region_df, region_sample_df):
    for _ in region_df.index:
        samples = list(region_sample_df.loc[_]['SAMPLE'])
        if len(samples) > 1:
            region_df.loc[_, 'patients'] = ' '.join(samples)
        elif len(samples) == 1:
            region_df.loc[_, 'patients'] = '{}'.format(samples[0])
    return region_df


def expand_samples(df):
    func = lambda x: pd.Series([i for i in x.split()])
    split = df['patients'].apply(func)
    df.drop('patients', axis=1, inplace=True)
    df = df.join(split.fillna(''))
    return df


def parse_sangers_for_seriereport(serie, db=None):
    all_sangers = get_sangers_serie_from_db(serie, db)
    all_sangers_df, failed_samples = combine_all_sangerdb_to_df(all_sangers)
    if not all_sangers_df.empty:
        all_sangers_df = get_sangerfrag_min_max(all_sangers_df)
        dfregions = all_sangers_df[['min', 'max']].astype(int).copy().drop_duplicates()
        dfregions['size'] = dfregions['max'] - dfregions['min']
        dfregions['patients'] = ''    
        df = group_samples_by_sanger(dfregions, all_sangers_df)
        df = expand_samples(df)
        df = df.join(all_sangers_df['gen']).drop_duplicates()
        return df, failed_samples
    elif all_sangers_df.empty:
        return pd.DataFrame(), failed_samples


def snpcheck2db(sample, serie, data, db):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = '''INSERT INTO snpcheck
    (SAMPLE, SERIE, DATA)
    VALUES ('{}', '{}', '{}')
    '''.format(sample, serie, data)
    c.execute(sql)
    conn.commit()
    conn.close()
    return


def get_snpcheck_serie(serie, db):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    sql = '''SELECT SAMPLE, DATA
    FROM snpcheck
    WHERE SERIE='{}'
    '''.format(serie)
    c.execute(sql)
    db_out = c.fetchall()
    conn.close()
    snpcheck_serie = dict()
    for _ in db_out:
        sample, data = _
        data = json.loads(data)
        snpcheck_serie[sample] = data
    return snpcheck_serie


def compare_snpchecks_serie(snpcheck_serie):
    dupcheck = defaultdict(list)
    for sample, data in snpcheck_serie.items():
        calls = list()
        for locus, _result in sorted(data['COMP'].items()):
            calls.append(data['NGS'][locus])
        calls = ':'.join(calls)
        dupcheck[calls].append(sample)
    return dupcheck


def mean_std_2db(mean, std, serie, sample, target, db):
    sql = '''INSERT INTO mean
    VALUES ('{}', '{}', '{}', '{}','{}')
    '''.format(sample, target, serie, mean, std)
    conn = sqlite3.connect(db)
    c = conn.cursor()
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        print(e)
    else:
        conn.commit()
    conn.close()

def riskscore_and_genotypes_2db(score, genotypes, serie, sample, target, db):
    sql = '''INSERT INTO riskscore
    VALUES ('{}', '{}', '{}', '{}', '{}')
    '''.format(sample, serie, target, score, genotypes)
    conn = sqlite3.connect(db)
    c = conn.cursor()
    try:
        c.execute(sql)
    except sqlite3.IntegrityError as e:
        print(e)
    else:
        conn.commit()
    conn.close()    