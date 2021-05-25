#!/bin/bash

mkdir -p vkr/build
cd vkr/build
cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=make ..
make