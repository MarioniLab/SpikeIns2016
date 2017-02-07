set -e
set -u

for f in Calero/trial_20160113 Calero/trial_20160325 Liora/test_20160906 Liora/test_20170201
do
    cd $f/analysis
    echo "knitr::knit('techanal.Rmd')" | R --no-save --vanilla;
    cd -
done
