echo "Setting up analysis environment"
wget https://raw.githubusercontent.com/banijolly/Genepi/master/covid19-environment.yml
conda env create -f covid19-environment.yml

source ~/anaconda3/etc/profile.d/conda.sh 
conda activate covid19-genepi

git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
pip install .

cd ..

echo "Enviroment set up complete"
