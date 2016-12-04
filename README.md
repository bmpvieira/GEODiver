# GeoDiver

[![Build Status](https://travis-ci.org/GeoDiver/GEODiver.svg?branch=master)](https://travis-ci.org/GeoDiver/GEODiver)
[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/GeoDiver/GEODiver/badges/quality-score.png?b=master)](https://scrutinizer-ci.com/g/GeoDiver/GEODiver/?branch=master)




## Introduction

GeoDiver is a web app that allows users to easily analyse GEO datasets.







## Installation
### Installation Requirements
* Ruby (>= 2.0.0)
  * Recommended to use rvm to install ruby
* R (=3.2.2)
  * Recommended to use R to install R

###Google API Setup
You will need need a GOOGLE_SECRET and GOOGLE_KEY in order to use the login system that GeoDiver utilises.

1. Go to 'https://console.developers.google.com'
2. Select your project.
3. Click 'Enable and manage APIs'.
4. Make sure "Contacts API" and "Google+ API" are on.
5. Go to Credentials, then select the "OAuth consent screen" tab on top, and provide an 'EMAIL ADDRESS' and a 'PRODUCT NAME'
6. Wait 10 minutes for changes to take effect.


### Installation
Simply run the following command in the terminal.

```bash
# Clone the repository.
git clone https://github.com/GeoDiver/GEODiver

# Move into GeoDiver source directory.
cd GEODiver

# Install R dependencies & Build and install the latest version of the webapp.
rake install 

# Start the web app
passenger start --envvar GOOGLE_KEY= --envvar GOOGLE_SECRET= -p 9292 -e production --sticky-sessions -d
```

##### Running From Source (Not Recommended)
It is also possible to run from source. However, this is not recommended.

```bash
# After cloning the web app and moving into the source directory 
# Install bundler
gem install bundler

# Use bundler to install dependencies
bundle install

# Optional: run tests and build the gem from source
bundle exec rake

# Run GeoDiver
bundle exec passenger start -h
# note that `bundle exec` executes GeoDiver in the context of the bundle

# Alternatively run Geodiver using the command line interface
bundle exec geodiver -h
```




## Launch GeoDiver

To configure and launch Geodiver, run the following from a command line from the GeoDiver root folder.

```bash
bundle exec passenger start -h

```
That's it! Open http://localhost:9292/ and start using GeoDiver!






## Advanced Usage

See `$ passenger start -h` for more information on all the options available when running GeoDiver.

# Config file
A Config file can be used to specify arguments - the default location of this file is in the home directory at `~/.geodiver.conf`. An examplar of the config file can be seen below.


```yaml
---
:num_threads: 8
:port: '9292'
:host: 0.0.0.0
:gd_public_dir: "/Users/ismailm/.geodiver"
:devel: true
```


<hr>

This program was developed at [QMUL](http://sbcs.qmul.ac.uk) as part of the [Bioinformatics Masters Course](http://www.qmul.ac.uk/postgraduate/taught/coursefinder/courses/121410.html).
