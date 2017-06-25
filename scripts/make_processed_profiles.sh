

REFERENCE=$1
SMOOTHSPAN=$2

echo "quanile-normalizing replication timing profiles"
quantile_normalize.sh 

echo "loess-smoothing replication timing profiles"
loess_smooth.R 300000
