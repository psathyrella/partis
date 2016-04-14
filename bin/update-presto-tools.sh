#!/bin/bash
# Script to update presto, changeo, alakazam, shazam and tigger
# 
# Author:  Jason Anthony Vander Heiden, Mohamed Uduman
# Date:    2016.04.13
# 
# Arguments:
#   -u = Bitbucket username. If not specified, username will be prompted.
#   -p = Bitbucket password. If not specified, password will be prompted.
#   -x = If specified, install the prototype topology, repertoire and presto
#        R packages.
#   -h = Display help.
#
# Creates Directories:
#   $outdir

# Set default parameters
PROTOTYPES=false
PASSWORD_PROMPT=true
USERNAME_PROMPT=true

outdir=$PWD/packages

# Print usage
usage () {
    echo "Usage: `basename $0` [OPTIONS]"
    echo "  -u  Bitbucket username."
    echo "  -p  Bitbucket password."
    echo "  -x  Install prototype packages."
    echo "  -h  This message."
}

# Get commandline arguments
while getopts "u:p:t:xh" OPT; do
    case "$OPT" in
    u)  USERNAME=$OPTARG
        USERNAME_PROMPT=false
        ;;
    p)  PASSWORD=$OPTARG
        PASSWORD_PROMPT=false
        ;;
    x)  PROTOTYPES=true
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# Ask for username if not provided
if $USERNAME_PROMPT; then
    read -p "Username: " USERNAME
fi

# Ask for password if not provided
if $PASSWORD_PROMPT; then
    read -s -p "Password: " PASSWORD
fi


# Create directories
mkdir -p $outdir

# Download bits
echo -e "\n\nDownloading default branch from repos..."
echo -e "================================================================================\n"
curl -# -f -u "${USERNAME}:${PASSWORD}" \
    https://bitbucket.org/kleinstein/presto/get/default.tar.gz -o $outdir/presto.tar.gz
if [ $? -ne 0 ]; then { echo -e "Download of presto failed." ; exit 1; } fi
curl -# -f -u "${USERNAME}:${PASSWORD}" \
    https://bitbucket.org/kleinstein/changeo/get/default.tar.gz -o $outdir/changeo.tar.gz
if [ $? -ne 0 ]; then { echo -e "Download of changeo failed." ; exit 1; } fi
curl -# -f -u "${USERNAME}:${PASSWORD}" \
    https://bitbucket.org/kleinstein/alakazam/get/default.tar.gz -o $outdir/alakazam.tar.gz
if [ $? -ne 0 ]; then { echo -e "Download of alakazam failed." ; exit 1; } fi
curl -# -f -u "${USERNAME}:${PASSWORD}" \
    https://bitbucket.org/kleinstein/shazam/get/default.tar.gz -o $outdir/shazam.tar.gz
if [ $? -ne 0 ]; then { echo -e "Download of shazam failed." ; exit 1; } fi
curl -# -f -u "${USERNAME}:${PASSWORD}" \
    https://bitbucket.org/kleinstein/tigger/get/default.tar.gz -o $outdir/tigger.tar.gz
if [ $? -ne 0 ]; then { echo -e "Download of tigger failed." ; exit 1; } fi

if $PROTOTYPES; then
    curl -# -f -u "${USERNAME}:${PASSWORD}" \
        https://bitbucket.org/javh/prototype-prestor/get/default.tar.gz \
        -o $outdir/prototype-prestor.tar.gz
    if [ $? -ne 0 ]; then { echo -e "Download of prototype-prestor failed." ; exit 1; } fi
    curl -# -f -u "${USERNAME}:${PASSWORD}" \
        https://bitbucket.org/javh/prototype-repertoire/get/default.tar.gz \
        -o $outdir/prototype-repertoire.tar.gz
    if [ $? -ne 0 ]; then { echo -e "Download of prototype-repertoire failed." ; exit 1; } fi
    curl -# -f -u "${USERNAME}:${PASSWORD}" \
        https://bitbucket.org/javh/prototype-topology/get/default.tar.gz \
        -o $outdir/prototype-topology.tar.gz
    if [ $? -ne 0 ]; then { echo -e "Download of prototype-topology failed." ; exit 1; } fi
fi


# Install python tools
echo -e "\n\nInstalling presto..."
echo -e "================================================================================\n"
mkdir -p $outdir/presto; cd $outdir/presto
tar -zxf $outdir/presto.tar.gz --wildcards --exclude="tests" --strip-components=1 \
    -C $outdir/presto
python3 setup.py install --user
cd -

echo -e "\n\nInstalling changeo..."
echo -e "================================================================================\n"
mkdir -p $outdir/changeo; cd $outdir/changeo
tar -zxf $outdir/changeo.tar.gz --wildcards --exclude="tests" --strip-components=1 \
    -C $outdir/changeo
python3 setup.py install --user
cd -


# Install R packages
RSCRIPT="options(repos=c(CRAN=\"http://watson.nci.nih.gov/cran_mirror\")); \
         library(roxygen2); library(devtools); \
         install_deps(\"${outdir}/package_directory\"); \
         document(\"${outdir}/package_directory\"); \
         build(\"${outdir}/package_directory\", vignettes=FALSE); \
         install(\"${outdir}/package_directory\")"

echo -e "\n\nInstalling alakazam..."
echo -e "================================================================================\n"
mkdir -p $outdir/alakazam
tar -zxf $outdir/alakazam.tar.gz --wildcards --exclude="tests" --strip-components=1 \
    -C $outdir/alakazam
Rscript -e "${RSCRIPT//package_directory/alakazam}"

echo -e "\n\nInstalling shazam..."
echo -e "================================================================================\n"
mkdir -p $outdir/shazam
tar -zxf $outdir/shazam.tar.gz --wildcards --exclude="tests" --strip-components=1 \
    -C $outdir/shazam
Rscript -e "${RSCRIPT//package_directory/shazam}"

echo -e "\n\nInstalling tigger..."
echo -e "================================================================================\n"
mkdir -p $outdir/tigger
tar -zxf $outdir/tigger.tar.gz --wildcards --exclude="tests" --strip-components=1 \
    -C $outdir/tigger
Rscript -e "${RSCRIPT//package_directory/tigger}"


# Install prototype R packages
if $PROTOTYPES; then
    echo -e "\n\nInstalling prototype-prestor..."
    echo -e "================================================================================\n"
    mkdir -p $outdir/prototype-prestor
    tar -zxf $outdir/prototype-prestor.tar.gz --wildcards --exclude="tests" \
        --strip-components=1 -C $outdir/prototype-prestor
    Rscript -e "${RSCRIPT//package_directory/prototype-prestor}"

    echo -e "\n\nInstalling prototype-repertoire..."
    echo -e "================================================================================\n"
    mkdir -p $outdir/prototype-repertoire
    tar -zxf $outdir/prototype-repertoire.tar.gz --wildcards --exclude="tests" \
        --strip-components=1 -C $outdir/prototype-repertoire
    Rscript -e "${RSCRIPT//package_directory/prototype-repertoire}"

    echo -e "\n\nInstalling prototype-topology..."
    echo -e "================================================================================\n"
    mkdir -p $outdir/prototype-topology
    tar -zxf $outdir/prototype-topology.tar.gz --wildcards --exclude="tests" \
        --strip-components=1 -C $outdir/prototype-topology
    Rscript -e "${RSCRIPT//package_directory/prototype-topology}"
fi
