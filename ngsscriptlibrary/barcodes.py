import os

def create_barcode_hussle_dict():
    """Turn TSV-file to dict.

    Read tab seperated file with our barcodenumber and the official one.
    Return dict with our number as key and the official as value.
    """
    basedir = os.path.dirname(os.path.abspath(__file__))
    barcodedict = dict()
    with open(os.path.join(basedir, 'docs/barcodehussle.txt')) as f:
        for line in f:
            if not line.split():
                continue
            else:
                us, them = line.strip().split()
                barcodedict[int(us)] = int(them)
    return barcodedict

def create_reverse_barcode_hussle_dict():
    """Turn TSV-file to dict.

    Read tab seperated file with our barcodenumber and the official one.
    Return dict with the official number as key and our number as value.
    """
    basedir = os.path.dirname(os.path.abspath(__file__))
    barcodedict = dict()
    with open(os.path.join(basedir, 'docs/barcodehussle.txt')) as f:
        for line in f:
            if not line.split():
                continue
            else:
                us, them = line.strip().split()
                barcodedict[int(them)] = int(us)
    return barcodedict    

def create_adapter_dict():
    """Turn TSV-file to dict.

    Read tab seperated file with barcode, seqI7 and seqI5.
    Return dict with barcode as key and a tuple of sequences as value.
    """
    basedir = os.path.dirname(os.path.abspath(__file__))
    adapterdict = dict()
    with open(os.path.join(basedir, 'docs/adapters.txt')) as f:
        for line in f:
            if not line.split():
                continue
            else:
                index, left, right = line.strip().split('\t')
                adapterdict[int(index)] = [left, right]
    return adapterdict

def create_reverse_adapter_dict():
    """Turn TSV-file to dict.

    Read tab seperated file with barcode, seqI7 and seqI5.
    Return dict with a tuple of sequences as key and a barcode as value.
    """
    basedir = os.path.dirname(os.path.abspath(__file__))
    adapterdict = dict()
    with open(os.path.join(basedir, 'docs/adapters.txt')) as f:
        for line in f:
            if not line.split():
                continue
            else:
                index, left, right = line.strip().split('\t')
                adapterdict[(left, right)] = int(index)
    return adapterdict

def create_adapter_code_dict():
    """Turn TSV-file to dict.

    Read tab seperated file with numbered seqI7 and seqI5.
    Return dict with sequence as key and number as value.
    """
    basedir = os.path.dirname(os.path.abspath(__file__))
    adapterdict = dict()
    with open(os.path.join(basedir, 'docs/adaptercodes.txt')) as f:
        for line in f:
            if not line.split():
                continue
            else:
                code, sequence = line.strip().split('\t')
                adapterdict[sequence] = code
    return adapterdict

