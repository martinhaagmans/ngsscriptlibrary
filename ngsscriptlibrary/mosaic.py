import os
import csv
import json
import glob
import sqlite3
from collections import namedtuple
from collections import defaultdict

import vcf
import pysam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class Mosaic:
    """Handle all things database for mosaic analysis."""

    def __init__(self, db):
        """Establisch connection with database and create a cursor."""

        self.db = db
        self.conn = sqlite3.connect(self.db)
        self.c = self.conn.cursor()

    def create_table(self):
        """Create a table for data per patient."""
        sql = '''CREATE TABLE IF NOT EXISTS patients
        (SAMPLE text NOT NULL,
        PAKKET text NOT NULL,
        ref text NOT NULL,
        var text NOT NULL,
        ref_below500 text NOT NULL,
        var_below500 text NOT NULL,
        PRIMARY KEY(SAMPLE, PAKKET))
        '''

        self.c.execute(sql)
        self.conn.commit()

    def create_reference_table(self, loci, reference):
        """Create a table with the refbase for all loci in the target."""
        ref = pd.DataFrame([(i,
                             pysam.faidx(reference,
                                         '{}:{}-{}'.format(i.split(':')[0],
                                                           i.split(':')[1],
                                                           i.split(':')[1])
                                         ).split('\n')[1])
                            for i in loci])
        ref.columns = ['Locus', 'ref']
        ref.set_index(['Locus'], inplace=True)
        ref.to_sql(con=self.conn, name='refbases', if_exists='replace')

    def add_data(self, data, sample, pakket):
        """Add dict with patient data to table."""
        ref = json.dumps(data['ref'])
        var = json.dumps(data['var'])
        ref_below500 = json.dumps(data['ref_below500'])
        var_below500 = json.dumps(data['var_below500'])
        sql = '''INSERT INTO patients
        VALUES ('{}', '{}', '{}', '{}', '{}', '{}')
        '''.format(sample, pakket, ref, var, ref_below500, var_below500)
        try:
            self.c.execute(sql)
        except sqlite3.IntegrityError as e:
            print(e)
        else:
            self.conn.commit()

    def get_reference_dict(self):
        """Query reftable for refbase per locus. Return a dict."""
        ref = dict()
        sql = 'SELECT * FROM refbases'
        self.c.execute(sql)
        for tup in self.c.fetchall():
            locus, base = tup
            ref[locus] = base
        return ref

    def get_archive_data(self, sample):
        """Query datatable for data for all samples. Return a dict."""
        sql = '''SELECT ref, var, ref_below500, var_below500
        FROM patients
        WHERE (SAMPLE !='{}')
        '''.format(sample)
        self.c.execute(sql)

        d = dict()
        d['allref'] = list()
        d['all_low_coverage'] = list()
        d['var'] = defaultdict(list)

        for _ in self.c.fetchall():
            ref, var, ref_below500, var_below500 = _

            ref = json.loads(ref)
            var = json.loads(var)

            ref_below500 = json.loads(ref_below500)
            var_below500 = json.loads(var_below500)

            [d['allref'].append(__) for __ in ref]
            [d['all_low_coverage'].append(__) for __ in ref_below500]
            [d['all_low_coverage'].append(__) for __ in var_below500.keys()]
            [d['var'][__].append(___) for __, ___ in var.items()]

        return d

    def get_sample_data(self, sample):
        """Query datatable for data for 1 sample. Return a dict."""
        d = dict()
        self.c.execute("SELECT var FROM patients WHERE (SAMPLE='{}')".format(sample))
        out = json.loads(self.c.fetchone()[0])
        for locus, data in out.items():
            d[locus] = dict()
            d[locus]['refperc'] = data[0]
            base_perc_list = data[1]
            for base_perc in base_perc_list:
                base, perc = base_perc
                d[locus][base] = perc
        return d

    def get_sample_low_coverage_var(self, sample):
        """Query datatable for low coverage variant bases for 1 sample. Return a dict."""
        self.c.execute("SELECT var_below500 FROM patients WHERE (SAMPLE='{}')".format(sample))
        return json.loads(self.c.fetchone()[0])

    def get_sample_low_coverage_ref(self, sample):
        """Query datatable for low coverage reference bases for 1 sample. Return a dict."""
        self.c.execute("SELECT ref_below500 FROM patients WHERE (SAMPLE='{}')".format(sample))
        return json.loads(self.c.fetchone()[0])

    def get_nrpatients(self):
        """Query datatable for nr of samples in table. Return a int."""
        self.c.execute("SELECT count(*) FROM patients")
        return int(self.c.fetchone()[0])

    @staticmethod
    def parse_sample_low_coverage_var(lc):
        lcd = dict()
        for locus in lc.keys():
            lcd[locus] = dict()
            base_perc_list = lc[locus][1]
            for base_perc in base_perc_list:
                base, perc = base_perc
                if base != 'I':
                    lcd[locus][base] = perc
                elif base == 'I':
                    lcd[locus][base] = perc[1]

        return lcd


def parse_bed(fn):
    with open(fn) as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            chrom, start, end, *_ = line
            yield chrom, start, end


def good_alignment(alignment):
    duplicate = alignment.is_duplicate
    mapped = not alignment.is_unmapped
    mapq20 = alignment.mapping_quality > 20
    return (not duplicate and mapped and mapq20)


def cigar_has_insertion(cigar):
    has_insertion = False
    if cigar is not None:
        if 'I' in cigar:
            has_insertion = True
    return has_insertion


def parse_cigartuple(cigartuple, read_start, chrom):
    for tup in cigartuple:
        operation, length = tup
        if operation == 2:
            continue
        elif operation == 1:
            locus = '{}:{}'.format(chrom, read_start)
            return locus, length
        else:
            read_start += length


def pos_in_interval(pos, intervalstart, intervalend):
    pos = int(pos)
    intervalstart = int(intervalstart)
    intervalend = int(intervalend)
    return pos >= intervalstart and pos <= intervalend


def get_indel_dicts(bamfile, target):
    samfile = pysam.AlignmentFile(bamfile, "rb")

    indel_coverage = defaultdict(int)
    indel_length = defaultdict(list)
    indel_length_coverage = dict()

    for c, s, e in parse_bed(target):
        s = int(s) - 151
        e = int(e) + 151

        for alignment in samfile.fetch(c, int(s), int(e)):

            if good_alignment(alignment) and cigar_has_insertion(alignment.cigarstring):

                read_start = alignment.get_reference_positions(full_length=True)[0]

                if read_start is None:
                    continue

                locus, length = parse_cigartuple(alignment.cigar, read_start,
                                                 alignment.reference_name)

                if pos_in_interval(locus.split(':')[1], s, e):

                    if locus in indel_length:
                        indel_length[locus].append(length)
                    else:
                        indel_length[locus] = [length]

                    indel_coverage[locus] += 1

    samfile.close()

    for locus, coverage in indel_coverage.items():
        indel_length_coverage[locus] = tuple(set(indel_length[locus])), int(coverage)

    return indel_length_coverage


def get_indel_dict_for_locus(bamfile, targetlocus):
    samfile = pysam.AlignmentFile(bamfile, "rb")

    indel_coverage = defaultdict(int)

    c, s = targetlocus.split(':')
    s = int(s) - 151
    e = int(s) + 151

    for alignment in samfile.fetch(c, int(s), int(e)):

        if good_alignment(alignment) and cigar_has_insertion(alignment.cigarstring):

            read_start = alignment.get_reference_positions(full_length=True)[0]

            if read_start is None:
                continue

            locus, _length = parse_cigartuple(alignment.cigar, read_start,
                                             alignment.reference_name)

            if str(targetlocus) == str(locus):
                indel_coverage[locus] += 1

    samfile.close()

    return indel_coverage


def parse_doc(fn, ref, loci):
    data = dict()
    with open(fn) as f:
        spamreader = csv.reader(f, delimiter='\t')
        _header = next(spamreader)
        for line in spamreader:
            nonref = list()
            locus, _TD, _ADS, DP, basecounts = line
            if locus not in loci:
                continue
            refbase = ref[locus]
            for bases in basecounts.split(' '):
                base, cov = bases.split(':')
                cov = int(cov)
                try:
                    basep = int(cov) / int(DP)
                except ZeroDivisionError:
                    basep = 0
                if refbase == base:
                    refp = basep
                elif refbase != base:
                    nonref.append((base, basep))
            locus_data = namedtuple('mosaic_out', 'DP, refpercentage, nonreflist')
            data[locus] = locus_data(int(DP), refp, nonref)

    return data


def parse_vcf(vcffile):
    vcfloci = list()
    for variant in vcf.Reader(open(vcffile), 'r'):
        if variant.FILTER:
            continue
        if variant.is_snp:
            locus = '{}:{}'.format(variant.CHROM, variant.POS)
        elif variant.is_deletion:
            locus = '{}:{}'.format(variant.CHROM, variant.POS + 1)
        elif variant.is_indel and not variant.is_deletion:
            locus = '{}:{}'.format(variant.CHROM, variant.POS)
        vcfloci.append(locus)
    return vcfloci


def split_data_for_database(data, vcfloci, DP_treshold=None):

    if DP_treshold is None:
        DP_treshold = 500

    d = dict()
    d['ref'] = list()
    d['var'] = dict()
    d['ref_below500'] = list()
    d['var_below500'] = dict()

    for l in data.keys():
        if l in vcfloci:
            continue
        if data[l].DP < DP_treshold:
            if data[l].refpercentage == 1:
                d['ref_below500'].append(l)
            elif data[l].refpercentage != 1:
                d['var_below500'][l] = data[l][1:3]

        elif data[l].DP >= DP_treshold:
            if data[l].refpercentage == 1:
                d['ref'].append(l)
            elif data[l].refpercentage != 1:
                d['var'][l] = data[l][1:3]

    return d


def get_mean(fn):
    df = pd.read_csv(fn, sep='\t')
    return df['Average_Depth_sample'].mean()


def bedfile_to_locilist(target):
    with open(target, 'r') as f:
        loci = list()
        for line in f:
            chrom, start, end, *_ = line.split()
            start = int(start)
            end = int(end)

            while start <= end:
                loci.append('{}:{}'.format(chrom, str(start)))
                start += 1

    return loci


def add_sampledata_to_database(bamfile, vcffile, docfile, sample, target, db):
    MDB = Mosaic(db)
    ref = MDB.get_reference_dict()
    loci = bedfile_to_locilist(target)
    sample_data = parse_doc(docfile, ref, loci)
    insertions = get_indel_dicts(bamfile, target)
    vcfloci = parse_vcf(vcffile)

    for l, data in sample_data.items():
        if l in insertions:
            sample_data[l].nonreflist.append(('I', (insertions[l][0],
                                                    insertions[l][1] / sample_data[l].DP)))
        elif l not in insertions:
            sample_data[l].nonreflist.append(('I', (0, 0)))

    data = split_data_for_database(sample_data, vcfloci)
    MDB.add_data(data, sample, 'SOv2')


def get_data_to_plot(sample, db):

    archive_plot_data = dict()
    sample_plot_data = dict()

    MDB = Mosaic(db)
    sample_data = MDB.get_sample_data(sample)
    archive_data = MDB.get_archive_data(sample)
    total_patients = MDB.get_nrpatients()

    for locus, data in sorted(sample_data.items()):

        ref_percentage = data['refperc']
        indel_percentage = data['I'][1]

        if ref_percentage < 0.99 or indel_percentage > 0.01:

            archive_plot_data[locus] = dict()
            sample_plot_data[locus] = data

            lc = archive_data['all_low_coverage'].count(locus)
            ref = archive_data['allref'].count(locus)
            loci_found = ref + lc

            if loci_found == total_patients - 1:
                for base in 'A C T G N D I'.split():
                    archive_plot_data[locus][base] = [0] * (total_patients - 1)

            elif loci_found < total_patients - 1:

                for base in 'A C T G N D I'.split():

                    if base not in archive_plot_data[locus]:
                        archive_plot_data[locus][base] = list()

                    for ref_base_perc in archive_data['var'][locus]:

                        for base_perc in ref_base_perc[-1]:

                            if base_perc[0] == base:

                                archive_plot_data[locus][base].append(base_perc[1])

            assert not loci_found > total_patients - 1

    return sample_plot_data, archive_plot_data


def plot_data(sample_plot_data, archive_plot_data, out, set_ylim=False):
    zoom = False
    sns.set_style('darkgrid')
    colordict = {'A': 'green', 'C': 'red', 'T': 'blue', 'G': 'black',
                 'D': 'sienna', 'I': 'orange'}

    fig = plt.figure(figsize=(12, 9))
    ax = plt.subplot()

    i = 0

    baselabels = list()
    boxplotdata = list()
    xticks = list()

    for locus, data in sorted(sample_plot_data.items()):

        for base, percentage in sorted(data.items()):
      
            if base == 'refperc' or base == 'N' or percentage == 0:
                continue
            if base == 'I':
                _length, percentage = percentage
                if percentage == 0:
                    continue

            if zoom is False and percentage > 0.1 and set_ylim is False:
                zoom = True
                    
            if base in baselabels:
                ax.plot(i, percentage, 'rD', c=colordict[base], markersize=10)
            elif base not in baselabels:
                ax.plot(i, percentage, 'rD', c=colordict[base], markersize=10, label=base)
                baselabels.append(base)

            if base != 'I':
                boxplotdata.append(archive_plot_data[locus][base])
            elif base == 'I':
                boxplotdata.append([_[1] for _ in archive_plot_data[locus][base]])

            xticks.append(locus)

            i += 1

    ax.boxplot(boxplotdata, notch=True, bootstrap=1000, meanline=True, showmeans=True,
               patch_artist=True, positions=[_ for _ in range(i)])

    for patch in ax.artists:
        r, g, b, _a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))

    plt.xticks([_ for _ in range(i)], xticks, rotation='vertical')
    plt.legend()
    if set_ylim:
        ax.set_ylim(0,0.1)
    fig.tight_layout()
    plt.savefig(out, dpi=120)
    return zoom


def parse_doc_for_literature_vars(docfile):
    data = dict()
    locus_data = namedtuple('mosaic_out', 'DP, A, C, G, T, D')
    with open(docfile) as f:
        spamreader = csv.reader(f, delimiter='\t')
        _header = next(spamreader)
        for line in spamreader:
            bp = dict()
            locus, _TD, _ADS, DP, basecounts = line

            for bases in basecounts.split(' '):
                base, cov = bases.split(':')
                cov = int(cov)
                try:
                    basep = int(cov) / int(DP)
                except ZeroDivisionError:
                    basep = 0

                bp[base] = basep

            data[locus] = locus_data(int(DP), bp['A'], bp['C'], bp['G'],
                                     bp['T'], bp['D'])

    return data


def get_known_mosaic_positions():
    litvars = dict()
    basedir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(basedir, 'docs', 'so_literature_vars.csv')) as f:
        fin = csv.reader(f)
        for line in fin:
            locus, varbase, cposvarstring = line
            litvars[locus] = varbase, cposvarstring
    return litvars