#!/bin/bash

# Defining colors
w=$(tput sgr0) 
r=$(tput setaf 1)
g=$(tput setaf 2) 
y=$(tput setaf 3) 
p=$(tput setaf 5) 
echo "${w}"
# Setting default values
kmer="" # --kmer-len in kraken2-build
min_l="" # --minimizer-len in kraken2-build
min_s="" # --minimizer-spaces in kraken2-build
read_len=75 # the read length in bracken-build
threads=`nproc`
#threads=$(($(nproc) - 2))
condapath=~/miniconda3
offline_files_folder=""
sudo_password=""

# Função auxiliar para rodar comandos com sudo usando expect
sudo_with_pass() {
    local cmd=$1
    expect <<EOF
        set timeout -1
        spawn bash -c "$cmd"
        expect {
            "*password*" {
                send "$sudo_password\r"
                exp_continue
            }
            eof
        }
EOF
}

#sudo_with_pass "sudo chown -R me:me /usr/local/lib/R/library"
# Function to display usage message
usage() {
    echo "Usage: $0 -p <condapath> -o <offline_files_folder> [-k <kmer>] [-m <minimizer_length>] [-s <minimizer_spaces>] [-r <read_length>] [-w <sudo_password>]"
    exit 1
}

# Checking if the required parameters are provided
if [ $# -lt 4 ]; then
    usage
fi

# Processing arguments
while getopts ":p:o:k:m:s:r:w:" option; do
    case "${option}" in
        p)
            condapath=${OPTARG}
            ;;
        o)
            offline_files_folder=${OPTARG}
            ;;
        k)
            kmer=${OPTARG}
            ;;
        m)
            min_l=${OPTARG}
            ;;
        s)
            min_s=${OPTARG}
            ;;
        r)
            read_len=${OPTARG}
            ;;
        w)
            sudo_password=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# Verifying if the required parameters are provided
if [ -z "$condapath" ] || [ -z "$offline_files_folder" ]; then
    usage
fi

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place
touch condaPath
echo "$condapath" > $dir/condaPath

source $condapath/etc/profile.d/conda.sh
sudo_with_pass "sudo apt-get update"
sudo_with_pass "sudo apt-get install libgeos-dev -y"
sudo_with_pass "sudo apt install libharfbuzz-dev libfribidi-dev libfreetype6-dev -y"
sudo_with_pass "sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev -y"
sudo_with_pass "sudo apt install libharfbuzz-dev rsync libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev pigz -y"
conda deactivate
echo 'installing conda environments...'
conda env create -f Installation/MTD.yml
conda env update -n MTD -f $dir/Installation/MTD_R_additions.yml
conda run -n MTD bash $dir/update_fix/Install.R.packages.MTD.sh
bash $dir/update_fix/check_R_pkg.MTD.sh

#$dir/update_fix/Install.R.packages.MTD.sh

sed -i 's/^rpy2[>=<]/# &/' $dir/Installation/pip.requirements
conda env create -f Installation/py2.yml
conda env create -f Installation/halla0820.yml
conda activate halla0820
#pip install --upgrade setuptools pip
#pip install -r Installation/pip.requirements
#pip install jenkspy matplotlib numpy pandas PyYAML scipy seaborn
#pip install --no-deps halla==0.8.20
conda deactivate
conda env create -f Installation/R412.yml
sed -i '/^# *rpy2/s/^# *//' $dir/Installation/pip.requirements
chmod +x $dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py
python3 -m py_compile $dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py

echo "${g}MTD installation progress:"
echo ">>                  [10%]${w}"

# conda activate py2 # install dependencies of py2 in case pip does work in conda yml
# pip install backports-functools-lru-cache==1.6.1 biom-format==2.0.1 cycler==0.10.0 h5py==2.10.0 hclust2==1.0.0 kiwisolver==1.1.0 matplotlib==2.2.5 numpy==1.16.6 pandas==0.24.2 pyparsing==2.4.7 pyqi==0.3.2 python-dateutil==2.8.1 pytz==2021.1 scipy==1.2.3 six==1.15.0 subprocess32==3.5.4
# conda deactivate

conda activate halla0820 # install dependencies of halla
#halla0820
conda install -n halla0820 -y -c conda-forge pkg-config
conda install -n halla0820 -y -c conda-forge ca-certificates openssl libcurl curl
conda install -n halla0820 -y -c conda-forge libuv
R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.22-7.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/mnormt_2.1.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/nlme_3.1-167.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/GPArotation_2024.3-1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/psych_2.5.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/foreign_0.8-89.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/R.methodsS3_1.8.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/R.oo_1.27.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/rtf_0.4-14.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('$dir/update_fix/pvr_pkg/psychTools_2.4.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('https://cran.r-project.org/src/contrib/XICOR_0.4.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('https://cran.r-project.org/src/contrib/mclust_6.1.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source')"
R -e "install.packages('$dir/update_fix/pvr_pkg/preprocessCore_1.72.0.tar.gz', repos=NULL, type='source')"
R -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
R -e 'remotes::install_url("https://cran.r-project.org/src/contrib/EnvStats_3.1.0.tar.gz", dependencies=TRUE)'
R -e 'remotes::install_version("Hmisc", version = "4.8-0", repos = "https://cloud.r-project.org")'
R -e "install.packages('https://cran.r-project.org/src/contrib/eva_0.2.6.tar.gz', repos=NULL, type='source')"

conda run -n halla0820 $dir/update_fix/check_R_pkg.halla0820.sh
conda deactivate
echo "${g}"
echo 'conda environments installed'
echo 'MTD installation progress:'
echo '>>>                 [15%]'
echo 'downloading virome database...'
echo "${w}"
conda activate MTD
sudo_with_pass "sudo apt-get update"
sudo_with_pass "sudo apt-get install rsync -y"
conda deactivate
##conda install -y python=3.10 
#conda install -n MTD -y -c bioconda metaphlan=4.0.6=pyhca03a8a_0 #Instalar no env MTD
conda activate MTD
conda install -n MTD -y -c conda-forge pkg-config
conda install -n MTD -y -c conda-forge ncbi-datasets-cli
conda install -n MTD -y -c conda-forge -c bioconda eggnog-mapper diamond

#Check if the file exists and have the same size before download
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/virushostdb.genomic.fna.gz
cp -f $offline_files_folder/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz .
#wget -c https://www.genome.jp/ftp/db/virushostdb/virushostdb.genomic.fna.gz

unpigz -f virushostdb.genomic.fna.gz
cat Installation/M33262_SIVMM239.fa virushostdb.genomic.fna > viruses4kraken.fa

echo "${g}"
echo "Adding additional viruses from NCBI Ref-Seq Viral to viruses4kraken.fa..."
echo "${w}"
cp $dir/manifest.virus.sh $offline_files_folder/Kraken2DB_micro/library/manifest.virus.sh
sed -i "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" $offline_files_folder/Kraken2DB_micro/library/manifest.virus.sh
$offline_files_folder/Kraken2DB_micro/library/manifest.virus.sh
sed -i -E 's/^>([^ ]+) (.+)/>\1 [\1] \2./' $offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna
#cat $offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna >> $dir/viruses4kraken.fa
# debug rsync error of kraken2-build
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/rsync_from_ncbi.pl
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
cp -f $dir/Installation/download_genomic_library.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh
cp -f $dir/Installation/download_genomic_library.sh $condapath/envs/MTD/libexec/download_genomic_library.sh
echo "${g}"
echo 'MTD installation progress:'
echo '>>>>                [20%]'
echo 'Preparing microbiome (virus, bacteria, archaea, protozoa, fungi, plasmid, UniVec_Core) database...'git 
echo "${w}"
# Kraken2 database building - Microbiome
#update for bacterial genomes
cp $dir/manifest.bacteria.sh $offline_files_folder/Kraken2DB_micro/library/manifest.bacteria.sh
sed -i "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" $offline_files_folder/Kraken2DB_micro/library/manifest.bacteria.sh
$offline_files_folder/Kraken2DB_micro/library/manifest.bacteria.sh

#Fix manifest.sh 
cp $dir/manifest.sh $offline_files_folder/Kraken2DB_micro/library/manifest.sh

# Substitui o path da pasta offline de instalacao no manifest.sh
sed -i "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" $offline_files_folder/Kraken2DB_micro/library/manifest.sh

DBNAME=kraken2DB_micro
echo "Downloading NCBI taxonomy database with Kraken2—please wait..."
# opção 1: tentar sem rsync
#kraken2-build --download-taxonomy --use-ftp --threads $threads --db "$DBNAME" $kmer $min_l $min_s
#kraken2-build --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
#Use a new script modified to optimze the download 
$dir/kraken2-build-download-taxonomy --download-taxonomy --threads $threads --db "$DBNAME" $kmer $min_l $min_s

# Use local files for archaea
echo "Downloading RefSeq Archaea library with Kraken2—please wait..."
#kraken2-build --use-ftp --download-library archaea --threads $threads --db $DBNAME $kmer $min_l $min_s
cp -f $dir/manifest.archea.sh $offline_files_folder/Kraken2DB_micro/library/manifest.archea.sh
sed -i "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" $offline_files_folder/Kraken2DB_micro/library/manifest.archea.sh
$offline_files_folder/Kraken2DB_micro/library/manifest.archea.sh

cp -f $dir/Installation/rsync_from_ncbi_archaea.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl

sed -i "13s|^.*|my \$local_download_dir = \"$offline_files_folder/Kraken2DB_micro/library/archaea/all/\";|" $condapath/envs/MTD/libexec/rsync_from_ncbi.pl

cp -f $offline_files_folder/Kraken2DB_micro/library/archaea/assembly_summary_archaea.txt $offline_files_folder/Kraken2DB_micro/library/archaea/assembly_summary.txt

chmod +x $condapath/envs/MTD/libexec/rsync_from_ncbi.pl

echo "Adding local archaeal sequences to Kraken2 database..."
kraken2-build --use-ftp --download-library archaea --threads $threads --db $DBNAME $kmer $min_l $min_s

# Restore original Kraken2 rsync script
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl

#Use local files for bacteria
cp -f $dir/Installation/rsync_from_ncbi_bacteria.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
sed -i "21s|^.*|my \$local_download_dir = \"$offline_files_folder/Kraken2DB_micro/library/bacteria/all/\";|" $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
cp -f $offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary_bacteria.txt $offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary.txt
chmod +x $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
echo "Adding local bacterial sequences to Kraken2 database..."
kraken2-build --use-ftp --download-library bacteria --threads $threads --db $DBNAME $kmer $min_l $min_s
cp -f $dir/Installation/rsync_from_ncbi.pl $condapath/envs/MTD/libexec/rsync_from_ncbi.pl
echo "Downloading RefSeq Protozoa library with Kraken2—please wait..."
kraken2-build --use-ftp --download-library protozoa --threads $threads --db $DBNAME $kmer $min_l $min_s
echo "Downloading RefSeq Fungi library with Kraken2—please wait..."
kraken2-build --use-ftp --download-library fungi --threads $threads --db $DBNAME $kmer $min_l $min_s

#Use local files for plasmid

#FIRST UPDATE PLASMID FILES
cp $dir/manifest.plasmid.sh $offline_files_folder/Kraken2DB_micro/library/manifest.plasmid.sh 
sed -i "s|^LOCAL_DIR=.*|LOCAL_DIR=$offline_files_folder/Kraken2DB_micro/library/plasmid/|" $offline_files_folder/Kraken2DB_micro/library/manifest.plasmid.sh 
$offline_files_folder/Kraken2DB_micro/library/manifest.plasmid.sh

sed -i "67s|^.*|    local_download_dir=\"$offline_files_folder/Kraken2DB_micro/library/plasmid/\"|" $dir/Installation/download_genomic_library_plasmid.sh
cp -f $dir/Installation/download_genomic_library_plasmid.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh    
cp -f $dir/Installation/download_genomic_library_plasmid.sh $condapath/envs/MTD/libexec/download_genomic_library.sh
echo "Adding local plasmid sequences to Kraken2 database..."
kraken2-build --download-library plasmid --threads $threads --db $DBNAME $kmer $min_l $min_s
cp -f $dir/Installation/download_genomic_library.sh $condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec/download_genomic_library.sh
cp -f $dir/Installation/download_genomic_library.sh $condapath/envs/MTD/libexec/download_genomic_library.sh

echo "Downloading UniVec_Core library with Kraken2—please wait..."
kraken2-build --use-ftp --download-library UniVec_Core --threads "$threads" --db "$DBNAME" $kmer $min_l $min_s

echo "Adding custom viral sequences (viruses4kraken.fa) to Kraken2 database..."
kraken2-build --add-to-library viruses4kraken.fa --threads $threads --db $DBNAME $kmer $min_l $min_s

echo "Building final Kraken2 database—this may take a while..."
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>              [30%]'
echo 'Preparing host (human) database...'
echo "${w}"
# Kraken2 database building - Human
DBNAME=kraken2DB_human
kraken2-build --use-ftp --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --use-ftp --download-library human --threads "$threads" --db "$DBNAME" $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>             [35%]'
echo 'Preparing host (mouse) database...'
echo "${w}"
# Kraken2 database building - Mouse
DBNAME=kraken2DB_mice
mkdir -p $DBNAME
cd $DBNAME
#wget -T 300 -t 5 -N --no-if-modified-since https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/GCF_000001635.27_GRCm39_genomic.fna.gz .
cp $offline_files_folder/GCF_000001635.27_GRCm39_genomic.fna.gz .

unpigz GCF_000001635.27_GRCm39_genomic.fna.gz
mv GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic.fa
cd ..
kraken2-build --use-ftp --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --add-to-library $DBNAME/GCF_000001635.27_GRCm39_genomic.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>            [40%]'
echo 'Preparing host (rhesus monkey) database...'
echo "${w}"
# Kraken2 database building - Rhesus macaque
DBNAME=kraken2DB_rhesus
mkdir -p $DBNAME
cd $DBNAME
#wget -T 300 -t 5 -N --no-if-modified-since https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/GCF_003339765.1_Mmul_10_genomic.fna.gz .
cp $offline_files_folder/GCF_003339765.1_Mmul_10_genomic.fna.gz .
unpigz GCF_003339765.1_Mmul_10_genomic.fna.gz
mv GCF_003339765.1_Mmul_10_genomic.fna GCF_003339765.1_Mmul_10_genomic.fa
cd ..
kraken2-build --use-ftp --download-taxonomy --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --add-to-library $DBNAME/GCF_003339765.1_Mmul_10_genomic.fa --threads $threads --db $DBNAME $kmer $min_l $min_s
kraken2-build --build --threads $threads --db $DBNAME $kmer $min_l $min_s

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>           [45%]'
echo 'Bracken database building...'
echo "${w}"
# Bracken database building
if [[ $kmer == "" ]]; then
    bracken-build -d $dir/kraken2DB_micro -t $threads -l $read_len
else
    bracken-build -d $dir/kraken2DB_micro -t $threads -l $read_len -k $kmer
fi

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>>>         [55%]'
echo 'installing HUMAnN3 databases...'
echo "${w}"
# install HUMAnN3 databases
mkdir -p $dir/HUMAnN/ref_database/
cd $dir/HUMAnN/ref_database/
#Link 403 forbidden
#wget -c http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.v296_201901.tar.gz
#Link working but slow
#wget -T 300 -t 5 -N --no-if-modified-since http://cmprod1.cibio.unitn.it/databases/HUMAnN/full_chocophlan.v296_201901.tar.gz
#Temporary cp solution
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/full_chocophlan.v296_201901.tar.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/HUMAnN_updated/full_chocophlan.v201901_v31.tar.gz .
cp $offline_files_folder/HUMAnN/full_chocophlan.v201901_v31.tar.gz .
#Link 403 forbidden
#wget -c http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz
#Link working but slow
#wget -T 300 -t 5 -N --no-if-modified-since http://cmprod1.cibio.unitn.it/databases/HUMAnN/uniref90_annotated_v201901.tar.gz
#Temporary cp solution
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/uniref90_annotated_v201901.tar.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/HUMAnN_updated/uniref90_annotated_v201901b_full.tar.gz .
cp $offline_files_folder/HUMAnN/uniref90_annotated_v201901b_full.tar.gz .

#Link 403 forbidden
#wget -c http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_v201901.tar.gz
#Link working but slow
#wget -c http://cmprod1.cibio.unitn.it/databases/HUMAnN/full_mapping_v201901.tar.gz
#Link source forge patrick-douglas
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/HUMAnN/ref_database/full_mapping_v201901.tar.gz
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/full_mapping_v201901.tar.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/HUMAnN_updated/full_mapping_v201901b.tar.gz .
cp $offline_files_folder/HUMAnN/full_mapping_v201901b.tar.gz .
mkdir -p $dir/HUMAnN/ref_database/chocophlan
#tar xzvf full_chocophlan.v296_201901.tar.gz -C chocophlan/
tar xzvf full_chocophlan.v201901_v31.tar.gz -C chocophlan/
#mkdir -p $dir/HUMAnN/ref_database/full_UniRef90
mkdir -p $dir/HUMAnN/ref_database/uniref
#tar xzvf uniref90_annotated_v201901.tar.gz -C full_UniRef90/
tar xzvf uniref90_annotated_v201901b_full.tar.gz -C uniref/
mkdir -p $dir/HUMAnN/ref_database/utility_mapping
#tar xzvf full_mapping_v201901.tar.gz -C utility_mapping/full_mapping_v201901b.tar.gz
tar xzvf full_mapping_v201901b.tar.gz -C utility_mapping/
cd $dir

humann_config --update database_folders nucleotide $dir/HUMAnN/ref_database/chocophlan
humann_config --update database_folders protein $dir/HUMAnN/ref_database/uniref
humann_config --update database_folders utility_mapping $dir/HUMAnN/ref_database/utility_mapping

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>      [70%]'
echo 'Fetching host (default: rhesus, human, mouse) references from local storage...'
echo "Local folder: $offline_files_folder"
echo "${w}"
# install host references
# download host GTF
    # download rhesus macaque GTF
#    wget -c http://ftp.ensembl.org/pub/release-104/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz -P ref_rhesus
#    wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz -P ref_rhesus
mkdir -p ref_rhesus && cp $offline_files_folder/Ref_genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz ref_rhesus

    # download human GTF
#    wget -c http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz -P ref_human
#    wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_human/Homo_sapiens.GRCh38.104.gtf.gz -P ref_human
mkdir -p ref_human && cp $offline_files_folder/Ref_genomes/Homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz ref_human

    # download mouse GTF
#    wget -c http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz -P ref_mouse
#    wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz -P ref_mouse
mkdir -p ref_mouse && cp $offline_files_folder/Ref_genomes/Mus_musculus/Mus_musculus.GRCm39.104.gtf.gz ref_mouse

# Building indexes for hisat2
echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>     [75%]'
echo 'Building host indexes (rhesus monkey) for hisat2...'
echo "${w}"
# rhesus macaques
mkdir -p hisat2_index_rhesus
cd hisat2_index_rhesus
cp ../ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz .
gzip -d Macaca_mulatta.Mmul_10.104.gtf.gz
mv Macaca_mulatta.Mmul_10.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
#wget -c http://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz #use ensembl genome to compatible with featureCount
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_rhesus/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
cp $offline_files_folder/Ref_genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz .
gzip -d Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
mv Macaca_mulatta.Mmul_10.dna.toplevel.fa genome.fa
hisat2-build --large-index -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>    [80%]'
echo 'Building host indexes (mouse) for hisat2...'
echo "${w}"
# mouse
mkdir -p hisat2_index_mouse
cd hisat2_index_mouse
cp ../ref_mouse/Mus_musculus.GRCm39.104.gtf.gz .
gzip -d Mus_musculus.GRCm39.104.gtf.gz
mv Mus_musculus.GRCm39.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
#wget -c http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz #use ensembl genome to compatible with featureCount
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
cp $offline_files_folder/Ref_genomes/Mus_musculus/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz .
gzip -d Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna.primary_assembly.fa genome.fa
hisat2-build --large-index -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>>   [85%]'
echo 'Building host indexes (human) for hisat2...'
echo "${w}"
# human
mkdir -p hisat2_index_human
cd hisat2_index_human
cp ../ref_human/Homo_sapiens.GRCh38.104.gtf.gz .
gzip -d Homo_sapiens.GRCh38.104.gtf.gz
mv Homo_sapiens.GRCh38.104.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
#wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz #use ensembl genome to compatible with featureCount
#wget -T 300 -t 5 -N --no-if-modified-since https://master.dl.sourceforge.net/project/mtd/MTD/ref_human/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
cp $offline_files_folder/Ref_genomes/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz .
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
hisat2-build --large-index -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

# # download preduild index for hisat2 # prebuild is from NCBI, may be not compatiable with featureCount
#     # H. sapiens
#     mkdir -p hisat2_index_human
#     cd hisat2_index_human
#     wget https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz
#     pigz -dc grch38_tran.tar.gz | tar xf -
#     cd ..

echo "${g}"
echo "Create a BLAST database for Magic-BLAST"
echo "${w}"
makeblastdb -in $dir/hisat2_index_human/genome.fa -dbtype nucl -parse_seqids -out $dir/human_blastdb/human_blastdb
makeblastdb -in $dir/hisat2_index_mouse/genome.fa -dbtype nucl -parse_seqids -out $dir/mouse_blastdb/mouse_blastdb
makeblastdb -in $dir/hisat2_index_rhesus/genome.fa -dbtype nucl -parse_seqids -out $dir/rhesus_blastdb/rhesus_blastdb

echo "${g}"
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%]'
echo 'installing R packages...'
echo "${w}"

# install R packages
conda deactivate
conda install -n py2 -y -c conda-forge pkg-config
conda activate R412
#$dir/update_fix/update_conda_pkgs.sh
conda install -n R412 -y -c conda-forge pkg-config

#conda run -n R412 bash $dir/update_fix/Install.R.packages.R412.sh

# debug in case libcurl cannot be located in the conda R environment
wget -T 300 -t 5 -N --no-if-modified-since https://cran.r-project.org/src/contrib/Archive/curl/curl_4.3.2.tar.gz
# if /usr/lib/x86_64-linux-gnu/pkgconfig/libcurl.pc exists, use it
if [ -f /usr/lib/x86_64-linux-gnu/pkgconfig/libcurl.pc ]; then
    locate_lib=/usr/lib/x86_64-linux-gnu/pkgconfig
    else 
    locate_lib=$(dirname $(locate libcurl | grep '\.pc'))
fi

#Install R412 env packages
cd $dir
conda deactivate
conda env remove -n R412 -y
rm -rf ~/miniconda3/envs/R412
rm -rf ~/miniconda3/pkgs/r-base-4.1.2-hde4fec0_0
rm -f ~/miniconda3/pkgs/r-base-4.1.2-hde4fec0_0*.tar.bz2
rm -f ~/miniconda3/pkgs/r-base-4.1.2-hde4fec0_0*.conda
conda clean --packages --tarballs -y
conda config --set channel_priority strict
conda env create -f $dir/Installation/R412.yml
conda activate R412
$dir/update_fix/Install.R.packages.R412_optimized.sh
Rscript -e 'install.packages("https://bioconductor.org/packages/3.19/bioc/src/contrib/UCSC.utils_1.0.0.tar.gz", repos=NULL, type="source", dependencies=FALSE); library(UCSC.utils); packageVersion("UCSC.utils")'
$dir/update_fix/check_R_pkg.R412.sh
conda deactivate
#Install Annotation tools for base enviroment
R -e 'BiocManager::install("GenomeInfoDb")'
bash $dir/update_fix/Install.R.AnnotPackages.base.sh
echo "${g}"
echo "*********************************"
echo "R packages version for conda envs"
echo "*********************************"
echo "${g}"
conda run -n MTD $dir/update_fix/check_R_pkg.MTD.sh
conda run -n R412 $dir/update_fix/check_R_pkg.R412.sh
conda run -n halla0820 $dir/update_fix/check_R_pkg.halla0820.sh
echo "${g}"
echo "*********************************"
echo ""
echo 'MTD installation progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%]'
echo "MTD installation is finished"
echo "${w}"
