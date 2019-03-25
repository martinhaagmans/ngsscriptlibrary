"""Annotate with UCSC Genome Browser or other sql RefGene database."""

__author__ = "Martin Haagmans (https://github.com/zaag)"
__license__ = "MIT"

import os

import pymysql
import pybedtools


class Annotation:
    """Annotate gene or region with RefGene.

    Method to query genomic intervals returns list of genenames.
    Method to query gene names returns list of genomic intervals.
    """

    def __init__(self, host=None, user=None, passwd=None, db=None,
                 genefile=None, genes=None):
        """Establish database connection."""
        if host is None and user is None and passwd is None:
            host = 'genome-mysql.cse.ucsc.edu'
            user = 'genome'
        if passwd is None:
            passwd = ''
        if db is None:
            db = 'hg19'

        self.host = host
        self.user = user
        self.db = db
        self.passwd = passwd
        self.conn = pymysql.connect(host=self.host, user=self.user,
                                    passwd=self.passwd, db=self.db)
        self.c = self.conn.cursor()
        self.genes = None

    def __del__(self):
        """Close database connection."""
        self.conn.close()

    def _parse_ucsc_output(self, out):
        """Read database output into list and return list."""
        if len(out) == 1:
            return [out[0][0]]
        elif len(out) > 1:
            return [tup[0] for tup in out]
        elif len(out) == [0]:
            ['NOGENE']

    def get_genename(self, chromosome, start=None, end=None, tx=True, cds=False):
        """Read locus or interval and return gene(s).

        Format locus/interval can be chr1:123456/chr1:123456-123789 or
        chromosome, start and (optional) end.
        Return a list with any gene(s) in the interval.
        """
        if cds:
            tx = False
        if start is None and end is None:
            if ':' not in chromosome:
                raise IOError('Geen chromosoom en start opgegeven')
            chromosome, region = chromosome.split(':')
            if '-' in region:
                start, end = region.split('-')
            else:
                start = region
        basetxsql = '''SELECT DISTINCT name2
        FROM refGene
        WHERE ((chrom='{c}') AND
        '''.format(c=chromosome)
        basecdssql = '''SELECT DISTINCT name2
        FROM refGene
        WHERE ((chrom='{c}') AND
        '''.format(c=chromosome)
        if end is not None:
            txsql = '''(('{s}'<txEnd AND '{s}'>txStart) AND
            ('{e}'<txEnd AND '{e}'>txStart)))
            '''.format(s=start, e=end)
            cdssql = '''(('{s}'<cdsEnd AND '{s}'>cdsStart) AND
            ('{e}'<cdsEnd AND '{e}'>cdsStart)))
            '''.format(s=start, e=end)
            txsql = basetxsql + txsql
            cdssql = basecdssql + cdssql
        elif end is None:
            txsql = '''('{s}'<txEnd AND '{s}'>txStart))
            '''.format(s=start)
            cdssql = '''('{s}'<cdsEnd AND '{s}'>cdsStart))Locus can be
            '''.format(s=start)
            txsql = basetxsql + txsql
            cdssql = basecdssql + cdssql

        if tx and not cds:
            self.c.execute(txsql)
        elif cds and not tx:
            self.c.execute(cdssql)

        if self.c.rowcount != 0:
            return self._parse_ucsc_output(self.c.fetchall())
        elif self.c.rowcount == 0:
            if end is not None:
                txsql = '''(('{s}'<txEnd AND '{s}'>txStart) OR
                ('{e}'<txEnd AND '{e}'>txStart)))
                '''.format(s=start, e=end)
                cdssql = '''(('{s}'<cdsEnd AND '{s}'>cdsStart) OR
                ('{e}'<cdsEnd AND '{e}'>cdsStart)))
                '''.format(s=start, e=end)
                txsql = basetxsql + txsql
                cdssql = basecdssql + cdssql
                if tx and not cds:
                    self.c.execute(txsql)
                elif cds and not tx:
                    self.c.execute(cdssql)
                if self.c.rowcount == 0:
                    return ['NOGENE']
                else:
                    return self._parse_ucsc_output(self.c.fetchall())
            elif end is None:
                return 'NOGENE'

    def get_region(self, gene, tx=True, cds=False):
        """Read gene name2 and return list of tuples.

        Parse gene name2 and return list of tuples with coding (cds) or
        transcription (tx) region for all NM-numbers of that gene.
        """
        if cds:
            tx = False

        txsql = """SELECT DISTINCT chrom, txStart, txEnd
        FROM refGene
        WHERE name2='{g}'
        """.format(g=gene)

        cdssql = """SELECT DISTINCT chrom, cdsStart, cdsEnd
        FROM refGene
        WHERE name2='{g}'
        """.format(g=gene)

        if tx and not cds:
            self.c.execute(txsql)
        elif cds and not tx:
            self.c.execute(cdssql)
        generegions = list()
        for i in self.c.fetchall():
            region = i[0], i[1], i[2]
            generegions.append(region)
        return generegions


class TargetAnnotation(Annotation):
    """Read BED file and optional (file with) list of genes.

    Method to annotate BED with gene names via RefSeq database
    Method to compare requested gene(s) with gene(s) in BED
    """

    def __init__(self, bedfile, genes=None, skip=25, **kwds):
        """Establish connection with database and parse BED-file and genelist.

        Raise IOerror if genes is not a file or list.
        """
        super().__init__(**kwds)

        self.skip = skip
        self.bedlocation = bedfile
        self.is_annotated = self._is_annotated()
        self.bed = self._parse_bed()

        if genes is None:
            self.genes = list()
        elif genes is not None:
            if isinstance(genes, list):
                genes = set(genes)
                self.genes = [gene for gene in genes]
            elif os.path.isfile(genes):
                self.genes = self._parse_genefile(genes)
        else:
            raise IOError("Can't understand genes. Should be list or file")

    def _is_annotated(self):
        """Check if there are >3 columns in BED. Return boolean."""
        with open(self.bedlocation, 'r') as f:
            cols = len(next(f).split())
        if cols == 3:
            return False
        elif cols > 3:
            return True

    def _parse_bed(self):
        """Read BED. Return list of lists."""
        bed = list()

        with open(self.bedlocation, 'r') as f:
            lines = [line.rstrip().split() for line in f]
        lines = list(line for line in lines if line)
        lines = sorted(lines, key=lambda x: (x[0], int(x[1]), int(x[2])))

        if not self.is_annotated:
            for line in lines:
                chromosome, start, end = line
                bed.append([chromosome, start, end])

        elif self.is_annotated:
            for line in lines:
                chromosome, start, end, gene = line[:4]
                bed.append([chromosome, start, end, gene])

        return bed

    def _parse_genefile(self, gf):
        """Read genefile into list and return list."""
        with open(gf, 'r') as f:
            lines = [line.rstrip() for line in f]
        lines = list(line for line in lines if line)
        return [l.split()[0] for l in lines]

    def _parse_annot_out(self, gl):
        if len(gl) == 1:
            return gl
        elif len(gl) > 1:
            return [_ for _ in gl]

    def annotate_bed(self):
        """Query RefSeq for every target in BED. Return list of lists."""
        annotated_bed = list()
        for target in self.bed:
            chromosome, start, end = target[:3]
            genename = self.get_genename(chromosome, int(start) + self.skip, int(end) - self.skip)
            annotated_bed.append([chromosome, start, end, '/'.join(genename)])
        self.is_annotated = True
        return annotated_bed

    def annotate_bed_and_filter_genes(self):
        """Annotate targets in BED and filter output with genelist.

        Query RefSeq for every target in BED, filter results with list of genes
        if the query returns >1 gene. Return list of lists.
        """
        annotated_bed = list()
        for target in self.bed:
            chromosome, start, end = target[:3]
            genename = self.get_genename(chromosome, int(start) + self.skip, int(end) - self.skip)
            if len(genename) > 1:
                _ = [gn for gn in genename if gn in self.genes]
                if len(_) >= 1:
                    genename = _
            annotated_bed.append([chromosome, start, end, '/'.join(genename)])
        self.is_annotated = True
        return annotated_bed

    def get_genes_in_bed(self):
        """Query RefSeq for every target in non-annotated BED into list and return list."""
        genesout = set()
        if not self.is_annotated:
            for target in self.bed:
                chromosome, start, end = target
                genename = self.get_genename(chromosome, int(start) + self.skip, int(end) - self.skip)
                genesout.update(self._parse_annot_out(genename))
            return list(genesout)
        elif self.is_annotated:
            return self.get_genes_in_annotated_bed(self.bed)

    def get_genes_in_annotated_bed(self, annotated_bed):
        """Collect genes from annotated BED into list and return list."""
        genesout = set()
        for target in annotated_bed:
            _chromosome, _start, _end, gene = target
            genesout.update([i for i in gene.split('/')])
        return list(genesout)

    def get_genes_not_found(self, bedgenes=None):
        """Compare genes in BED with genes requested  into list and return list."""
        if bedgenes is None:
            bedgenes = self.get_genes_in_bed()
        notfound = [gene for gene in self.genes if gene not in bedgenes]
        if len(notfound) > 0:
            notfound.sort()
        return notfound

    def get_genes_not_requested(self, genesout=None):
        """Compare genes requested with genes in BED into list and return list."""
        if genesout is None:
            genesout = self.get_genes_in_bed()
        notrequested = [gene for gene in genesout if gene not in self.genes]
        return notrequested

    def report_genecomp(self):
        """Report differences between genes in BED and genes requested."""
        if not self.is_annotated:
            annotated_bed = self.annotate_bed()
        elif self.is_annotated:
            annotated_bed = self.bed
        genes_in_bed = self.get_genes_in_annotated_bed(annotated_bed)
        genes_not_found = self.get_genes_not_found(bedgenes=genes_in_bed)
        genes_not_requested = self.get_genes_not_requested(genesout=genes_in_bed)
        return genes_not_found, genes_not_requested


def annotate_bed(bedfile, output, genes=None, merge=True):
    """Merge and annotate bedfile.

    A file with genenames is optional and used to filter if needed.
    Write merged and annotated BED-file to output.
    """
    if merge:
        bed = pybedtools.BedTool(bedfile)
        bed = bed.sort()
        bed_merged = bed.merge()
        bed_merged.saveas(bedfile)

    if genes:
        TA = TargetAnnotation(bedfile, genes=genes)
        bed_annotated = TA.annotate_bed_and_filter_genes()
    elif genes is None:
                TA = TargetAnnotation(bedfile)
                bed_annotated = TA.annotate_bed()

    with open(output, 'w') as f:
        for line in bed_annotated:
            chromosome, start, end, gene = line
            f.write('{}\t{}\t{}\t{}\n'.format(chromosome, start, end, gene))
