#!/usr/local/bin/bash

# source the correct python
source /o/ishivvers/my_python/bin/activate

# start the mongodb
numactl --interleave=all mongod --config /o/ishivvers/db/mongodb.conf &

# start a bunch of flask instances
for port in `seq 0 3`
  do
  python /o/ishivvers/zeropoint_server/run_photozpe.py $port &
  done
echo servers started successfully.
