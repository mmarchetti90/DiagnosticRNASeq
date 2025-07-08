#!/bin/bash

dockerfile=$1
name=$2
tag=$3

mkdir -p sif_images
mkdir -p tar_images
docker build -t ${name}:${tag} -f ${dockerfile} . &> docker_${name}_${tag}.log
docker save -o ${name}_${tag}.tar localhost/${name}:${tag}
singularity build ${name}_${tag}.sif docker-archive://${name}_${tag}.tar
mv ${name}_${tag}.tar tar_images/
mv ${name}_${tag}.sif sif_images/
