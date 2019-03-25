# NGSSCRIPTLIBRARY

Collection of scripts to help with NGS data-analysis in the AMC.

## Prerequisites

zlib and BEDtools:

```bash
sudo apt-get install zlib-dev
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.0/bedtools-2.27.0.tar.gz
tar -zxvf bedtools-2.27.0.tar.gz
cd bedtools2
make
sudo cp bedtools2/bin/* /usr/local/bin
```

python3-tkinter:
```bash
sudo apt-get install python3-tk
```


## Install

```bash
sudo pip install ngsscriplibrary
```
