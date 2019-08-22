import os
import json
import sqlite3
import urllib.request


def download_file(url, file_name):
    urllib.request.urlretrieve(url, file_name)


def get_picard_header():
    picardheader = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'docs', 'picardheader.txt')
    out = list()
    with open(picardheader, 'r') as header:
        for line in header:
            out.append(line)
    return out


def boolean_to_number(x):
    if x == 'False':
        x = 0
    elif not x:
        x = 0
    elif x == 'True':
        x = 1
    elif x:
        x = 1
    return x


class TargetDatabase:
    """Establish connection with AMC's db of diagnostic tests.
    Download db from GitHub if location is not given.

    Raises OSError if given location is not a file.
    """
    def __init__(self, db=None):
        if db is None:
            url = ('https://github.com/martinhaagmans/ngstargets/'
                   'blob/master/varia/captures.sqlite?raw=true')
            HOME = os.path.expanduser('~')
            download_file(url, '{}/captures.sqlite'.format(HOME))
            self.db = '{}/captures.sqlite'.format(HOME)

        elif db is not None:
            self.db = db

        if not os.path.isfile(self.db):
            dbloc = os.path.abspath(self.db)
            raise OSError('Database file {} bestaat niet'.format(dbloc))

        self.conn = sqlite3.connect(self.db)
        self.c = self.conn.cursor()
        self.tables = ['genesis', 'aandoeningen', 'captures',
                       'pakketten', 'panels']

    def __repr__(self):
        return "{}(db='{}')".format(self.__class__.__name__,
                                    os.path.abspath(self.db))

    def __del__(self):
        self.conn.close()
        if self.db == '{}/captures.sqlite'.format(os.path.expanduser('~')):
            os.remove(self.db)

    @staticmethod
    def split_version(vname):
        return vname.split('v')

    @staticmethod
    def number_to_boolean(x):
        "Turn 0 to False, 1 to True. Return Boolean."
        if x == 0:
            return False
        elif x == 1:
            return True

    def list_db_output(self):
        """Create a list of sqlite output tuples.
        Return list.
        """
        return [val for tup in self.c.fetchall() for val in tup]

    def combine_version(self, name, versies):
        """Combine capture with highest version number from list.
        Return string.
        """
        versies.sort()
        hoogste = versies[-1]
        return '{}v{}'.format(name, hoogste)

    def get_current_capture(self, genesis):
        """Find and combine capture and version for given genesiscode.
        Return string.
        """
        sql = """SELECT capture
        FROM genesis
        WHERE genesis='{}'
        """.format(genesis)
        self.c.execute(sql)
        out = self.c.fetchone()
        sql = """SELECT versie
        FROM captures
        WHERE capture='{}'
        """.format(out[0])
        self.c.execute(sql)
        return self.combine_version(out[0], self.list_db_output())

    def get_current_pakket(self, genesis):
        """Find and combine pakket and version for given genesiscode.
        Return string.
        """
        sql = """SELECT pakket
        FROM genesis
        WHERE genesis='{}'
        """.format(genesis)
        self.c.execute(sql)
        out = self.c.fetchone()
        sql = """SELECT versie
        FROM pakketten
        WHERE pakket='{}'
        """.format(out[0])
        self.c.execute(sql)
        return self.combine_version(out[0], self.list_db_output())

    def get_current_panel(self, genesis):
        """Find and combine panel and version for given genesiscode.
        Return string or NoneType.
        """
        sql = """SELECT panel
        FROM genesis
        WHERE genesis='{}'
        """.format(genesis)
        self.c.execute(sql)
        out = self.c.fetchone()
        if out[0] == 'None':
            return None
        sql = """SELECT versie
        FROM panels
        WHERE panel='{}'
        """.format(out[0])
        self.c.execute(sql)
        panel = self.combine_version(out[0], self.list_db_output())
        return panel.replace('v', 'typeAv')

    def get_aandoening(self, genesis):
        sql = """SELECT aandoening
        FROM aandoeningen
        WHERE genesis='{}'
        """.format(genesis)
        self.c.execute(sql)
        return self.c.fetchone()[0]

    def do_x_in_pipeline(self, genesis, x):
        "Get 0/1 for x. Return boolean."
        sql = """SELECT {}
        FROM genesis
        WHERE genesis='{}'
        """.format(x, genesis)
        self.c.execute(sql)
        out = self.c.fetchone()
        return TargetDatabase.number_to_boolean(out[0])

    def check_if_genesis_exists(self, genesis):
        """Get all genesis codes from table and check if genesis is in list.
        Return boolean.
        """
        self.c.execute('SELECT genesis FROM genesis')
        return genesis in self.list_db_output()

    def get_info_for_genesis(self, genesis):
        info = dict()
        info['capture'] = self.get_current_capture(genesis)
        info['pakket'] = self.get_current_pakket(genesis)
        info['panel'] = self.get_current_panel(genesis)
        info['aandoening'] = self.get_aandoening(genesis)
        analyses = ['cnvscreening', 'cnvdiagnostiek', 'mozaiekdiagnostiek']
        for analysis in analyses:
            info[analysis] = self.do_x_in_pipeline(genesis, analysis)
        return info

    def get_todo_for_amplicon_confirmation(self, genesis):
            genesis = genesis.replace(' ', '')
            genesis = genesis.lower()
            todo_list = {'target': genesis, 
                         'amplicon': True, 
                         'capture': 'AMPLICON', 
                         'pakket': 'AMPLICON',
                         'panel': None,
                         'cnvscreening': False,
                         'cnvdiagnostiek': False,
                         'mozaiek': False,
                         'riskscore': False,
                         'capispakket': True,
                         }
            return todo_list

    def get_todo(self, genesis):
        """Get current capture, pakket, panel, and todo cnv and mosaic.
        Return dict.
        """
        if 'chr' in genesis.lower() and ':' in genesis:
            return self.get_todo_for_amplicon_confirmation(genesis)

        todo_list = dict()

        if not self.check_if_genesis_exists(genesis):
            raise IOError('{} does not exist in {}'.format(genesis, self.db))
        todo_list['amplicon'] = self.do_x_in_pipeline(genesis, 'amplicon')
        if todo_list['amplicon']:
            todo_list['capture'] = 'AMPLICON'
            todo_list['pakket'] = 'AMPLICON'
            todo_list['panel'] = None
        else:
            todo_list['capture'] = self.get_current_capture(genesis)
            todo_list['pakket'] = self.get_current_pakket(genesis)
            todo_list['panel'] = self.get_current_panel(genesis)
        todo_list['cnvscreening'] = self.do_x_in_pipeline(genesis, 'cnvscreening')
        todo_list['cnvdiagnostiek'] = self.do_x_in_pipeline(genesis,'cnvdiagnostiek')
        todo_list['mozaiek'] = self.do_x_in_pipeline(genesis,'mozaiekdiagnostiek')
        todo_list['riskscore'] = self.do_x_in_pipeline(genesis, 'riskscore')    
                                                                                                              
        if todo_list['capture'].split('v')[0] == todo_list['pakket'].split('v')[0]:
            todo_list['capispakket'] = True
        elif todo_list['capture'].split('v')[0] != todo_list['pakket'].split('v')[0]:
            todo_list['capispakket'] = False
            todo_list['pakketgenen'] = self.get_genes_for_vpakket(todo_list['pakket'])
        return todo_list

    def get_all_tests(self):
        "Return list of all available AMC diagnostic tests from database."
        self.c.execute('SELECT genesis FROM genesis')
        return self.list_db_output()

    def get_all_captures(self):
        "Retrieve all existing captures. Return list."
        sql = "SELECT DISTINCT capture FROM captures"
        self.c.execute(sql)
        return self.list_db_output()

    def get_all_pakketten(self):
        "Retrieve all existing pakketten. Return list."
        sql = "SELECT DISTINCT pakket FROM pakketten"
        self.c.execute(sql)
        return self.list_db_output()

    def get_all_panels(self):
        "Retrieve all existing panels. Return list."
        sql = "SELECT DISTINCT panel FROM panels"
        self.c.execute(sql)
        return self.list_db_output()

    def get_all_captures_for_genesis(self, genesis):
        "Retrieve existing capture(s) for test. Return list."
        self.check_if_genesis_exists(genesis)
        sql = '''SELECT DISTINCT capture
        FROM genesis WHERE (genesis='{}')
        '''.format(genesis)
        self.c.execute(sql)
        return self.list_db_output()

    def get_all_pakketten_for_genesis(self, genesis):
        "Retrieve existing capture(s) for test. Return list."
        self.check_if_genesis_exists(genesis)
        sql = '''SELECT DISTINCT pakket
        FROM genesis WHERE (genesis='{}')
        '''.format(genesis)
        self.c.execute(sql)
        return self.list_db_output()

    def get_all_panels_for_genesis(self, genesis):
        "Retrieve existing capture(s) for test. Return list."
        self.check_if_genesis_exists(genesis)
        sql = '''SELECT DISTINCT panel
        FROM genesis WHERE (genesis='{}')
        '''.format(genesis)
        self.c.execute(sql)
        return self.list_db_output()

    def get_aandoening_for_genesis(self, genesis):
        "Retrieve aandoening. Return string."
        self.check_if_genesis_exists(genesis)
        sql = '''SELECT DISTINCT aandoening
        FROM aandoeningen WHERE (genesis='{}')
        '''.format(genesis)
        self.c.execute(sql)
        return self.c.fetchone()[0]

    def get_all_versions_for_capture(self, capture):
        sql = '''SELECT DISTINCT versie
        FROM captures
        WHERE capture='{}'
        '''.format(capture)
        self.c.execute(sql)
        out = list()
        for versie in self.list_db_output():
            out.append('{}v{}'.format(capture, versie))
        return out

    def get_all_versions_for_pakket(self, pakket):
        sql = '''SELECT DISTINCT versie
        FROM pakketten
        WHERE pakket='{}'
        '''.format(pakket)
        self.c.execute(sql)
        out = list()
        for versie in self.list_db_output():
            out.append('{}v{}'.format(pakket, versie))
        return out

    def get_all_versions_for_panel(self, panel):
        sql = '''SELECT DISTINCT versie
        FROM panels
        WHERE panel='{}'
        '''.format(panel)
        self.c.execute(sql)
        out = list()
        for versie in self.list_db_output():
            out.append('{}typeAv{}'.format(panel, versie))
        return out

    def get_all_genesiscodes_for_vcapture(self, vcapture):
        """Retrieve genesis for vcapture. Return list."""
        capture, _version = self.split_version(vcapture)
        sql = '''SELECT genesis
        FROM genesis
        WHERE capture='{}'
        '''.format(capture)
        self.c.execute(sql)
        return self.list_db_output()

    def get_oid_for_vcapture(self, vcapture):
        """Retrieve OID. Return string.
        """
        capture, version = self.split_version(vcapture)
        sql = '''SELECT oid
        FROM captures
        WHERE (capture='{}' and versie={})
        '''.format(capture, int(version))
        self.c.execute(sql)
        oid = self.c.fetchone()[0]
        return str(oid)

    def get_all_info_for_vcapture(self, vcapture):
        """Retrieve OID, verdund and size.
        Return list of tuples
        """
        capture, version = self.split_version(vcapture)
        sql = '''SELECT oid, verdund, grootte
        FROM captures
        WHERE (capture='{}' and versie={})
        '''.format(capture, int(version))
        self.c.execute(sql)
        return self.c.fetchall()

    def get_genes_for_vcapture(self, vcapture):
        """Retrieve genes. Return list."""
        capture, version = self.split_version(vcapture)
        sql = '''SELECT genen
        FROM captures
        WHERE (capture='{}' and versie={})
        '''.format(capture, int(version))
        self.c.execute(sql)
        genes = self.c.fetchone()[0]
        return json.loads(genes)

    def get_genes_for_vpakket(self, vpakket):
        """Retrieve genes. Return list."""
        pakket, version = self.split_version(vpakket)
        sql = '''SELECT genen
        FROM pakketten
        WHERE (pakket='{}' and versie={})
        '''.format(pakket, int(version))
        self.c.execute(sql)
        genes = self.c.fetchone()[0]
        return json.loads(genes)

    def get_genes_for_vpanel(self, vpanel):
        """Retrieve genes. Return list."""
        panel, version = self.split_version(vpanel)
        sql = '''SELECT genen
        FROM panels
        WHERE (panel='{}' and versie={})
        '''.format(panel, int(version))
        self.c.execute(sql)
        genes = self.c.fetchone()[0]
        return json.loads(genes)

    def get_all_info_for_vpakket(self, vpakket):
        """Retrieve size and genes.
        Return list of tuples.
        """
        pakket, version = self.split_version(vpakket)
        sql = '''SELECT grootte, genen
        FROM pakketten
        WHERE (pakket='{}' and versie={})
        '''.format(pakket, int(version))
        self.c.execute(sql)
        dbout = self.c.fetchall()
        out = list()
        for pakket in dbout:
            grootte, genen, *_ = pakket
            genen = sorted(json.loads(genen))
            out.append((grootte, genen))
        return out

    def get_all_info_for_vpanel(self, vpanel):
        """Retrieve size and genes.
        Return list of tuples
        """
        if 'typeA' in vpanel:
            vpanel = vpanel.replace('typeA', '')
        panel, version = self.split_version(vpanel)
        sql = '''SELECT grootte, genen
        FROM panels
        WHERE (panel='{}' and versie={})
        '''.format(panel, int(version))
        self.c.execute(sql)
        dbout = self.c.fetchall()
        out = list()
        for panel in dbout:
            grootte, genen, *_ = panel
            genen = sorted(json.loads(genen))
            out.append((grootte, genen))
        return out

    def get_capture_for_pakket(self, pakket):
        """Retrieve capture and version.
        Return string.
        """
        sql = '''SELECT capture
        FROM genesis
        WHERE pakket='{}'
        '''.format(pakket)
        self.c.execute(sql)
        cap = self.c.fetchone()
        sql = """SELECT versie
        FROM captures
        WHERE capture='{}'
        """.format(cap[0])
        self.c.execute(sql)
        return self.combine_version(cap[0], self.list_db_output())

    def get_pakket_for_panel(self, panel):
        """Retrieve pakket and version.
        Return string.
        """
        sql = '''SELECT pakket
        FROM genesis
        WHERE panel='{}'
        '''.format(panel)
        self.c.execute(sql)
        pakket = self.c.fetchone()
        sql = """SELECT versie
        FROM pakketten
        WHERE pakket='{}'
        """.format(pakket[0])
        self.c.execute(sql)
        return self.combine_version(pakket[0], self.list_db_output())

    def get_amplicon_genesiscode(self):
        sql = """SELECT genesis
        FROM genesis
        WHERE amplicon=1 
        """
        self.c.execute(sql)
        return self.list_db_output()

    def change(self, sql):
        "Execute and commit sql statement."
        self.c.execute(sql)
        self.conn.commit()


class TargetFiles:

    def __init__(self, ngstargetdir):
        self.captures = os.path.join(ngstargetdir, 'captures')
        self.pakketten = os.path.join(ngstargetdir, 'pakketten')
        self.panels = os.path.join(ngstargetdir, 'panels')
        self.categories = ['capture', 'pakket', 'panel']

    def _check(self, category):
        if str(category) not in self.categories:
            raise IOError('Categorie moet capture, pakket of panel zijn')

    def _get_basedir(self, category):
        if category == 'capture':
            return self.captures
        elif category == 'pakket':
            return self.pakketten
        elif category == 'panel':
            return self.panels

    def get_sanger_target(self, name, category):
        self._check(category)
        fn = '{}_targets.bed'.format(name)
        sangertarget = os.path.join(self._get_basedir(category), fn)
        return sangertarget

    def get_variant_call_target(self, name, category):
        self._check(category)
        fn = '{}_generegions.bed'.format(name)
        varcalltarget = os.path.join(self._get_basedir(category), fn)
        return varcalltarget

    def get_picard_target(self, name):
        fn = '{}_target.picard.interval_list'.format(name)
        picardtarget = os.path.join(self.captures, fn)
        return picardtarget

    def get_cnv_target(self, name):
        fn = '{}_target.bed'.format(name)
        cnvtarget = os.path.join(self.captures, fn)
        return cnvtarget
