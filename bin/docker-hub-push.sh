#!/bin/bash

sudo docker build --tag psathyrella/partis .
# sudo docker tag <local tag> psathyrella/partis
sudo docker login -u psathyrella
sudo docker push psathyrella/partis
