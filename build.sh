#!/bin/bash
sudo rm -rf wheelhouse build pkgs rayx.cpp rayx/build
sudo docker run -it -v .:/io quay.io/pypa/manylinux_2_28_x86_64 /io/builder.sh