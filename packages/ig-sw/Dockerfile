#Dockerfile to build ig-sw container images
#Base image is ubuntu

FROM ubuntu:14.04

RUN apt-get update
RUN apt-get -y dist-upgrade

RUN apt-get install -y scons gcc g++ zlib1g-dev libpthread-stubs0-dev


ADD . /ig-sw/
WORKDIR "ig-sw/src/ig_align/"
RUN scons
WORKDIR "/ig-sw/src/"
RUN g++ test.cpp -lz
RUN ./a.out
WORKDIR "/ig-sw/"
