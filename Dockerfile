FROM ubuntu:14.04
MAINTAINER Connor McCoy

ADD ./debs /tmp/debs

RUN apt-get update -qq
RUN apt-get install -qq -y cython python-pip openjdk-7-jre-headless
RUN dpkg -i /tmp/debs/*.deb

CMD vdjalign --help && ighutil --help
