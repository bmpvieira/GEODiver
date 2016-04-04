# GEODiver
#
# VERSION 1.0.0

FROM ubuntu:16.04
MAINTAINER Bruno Vieira <mail@bmpvieira.com>

LABEL Description="GEODiver" Version="1.0.0"

ENV USER root
ENV TMPDIR /tmp

# Base packages
RUN apt-get update
RUN apt-get install -y build-essential autoconf automake curl wget git ssh vim npm python python-dev

# Install R and Ruby
RUN apt-get install -y ruby ruby-dev r-base r-base-dev

# Install other dependencies
RUN apt-get install -y libltdl-dev curl libcurl3 libxml2 libxml2-dev libcairo2 libcairo2-dev libxt-dev libxaw7 libxaw7-dev

# Add source code
ADD . /app
WORKDIR /app

# Build
RUN gem install bundler
RUN bundle install
RUN rake install

# Run
RUN passenger start --envvar GOOGLE_KEY= --envvar GOOGLE_SECRET= -p 9292 -e production --sticky-sessions -d
