#!/bin/bash

# Download and extract IgBLAST
VERSION="1.17.1"
ncbdir=packages/ncbi-igblast-${VERSION}
imdir=$ncbdir/immcantation/scripts
sh_igbdir=$ncbdir/igbl-db  # /~/share/igblast
sh_imgdir=$ncbdir/imgt  # ~/share/germlines/imgt
export PATH=$PATH:$PWD/$imdir

# wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${VERSION}/ncbi-igblast-${VERSION}-x64-linux.tar.gz
# tar -zxf ncbi-igblast-${VERSION}-x64-linux.tar.gz
# cp $ncbdir/bin/* ~/bin

# Download reference databases and setup IGDATA directory
$imdir/fetch_igblastdb.sh -o $sh_igbdir
cp -r $ncbdir/internal_data $sh_igbdir
cp -r $ncbdir/optional_file $sh_igbdir

# Build IgBLAST database from IMGT reference sequences
$imdir/fetch_imgtdb.sh -o $sh_imgdir
$imdir/imgt2igblast.sh -i $sh_imgdir -o $sh_igbdir
