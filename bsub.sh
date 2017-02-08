R CMD INSTALL --clean package

cd simulations/diffexp
bsub -R "rusage[mem=10000]" -n 10 -e "log.err" -o "log.out" R CMD BATCH --no-save desim.R
cd -

cd simulations/variance
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save varsim.R
cd -

cd simulations/clustering
bsub -R "rusage[mem=5000]" -n 1 -e "log.err" -o "log.out" R CMD BATCH --no-save clustsim.R
cd - 
