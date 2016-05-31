# ConcurrentDecoder [![Build Status](https://travis-ci.org/sifaserdarozen/ConcurrentDecoder.png)](https://travis-ci.org/sifaserdarozen/ConcurrentDecoder)
multi channel decoder experiment

### What is it?
Experiment code that uses Nvidia cuda in decoding g711 / g722 rtp payloads. 
v0.1 version corresponds to code of 2014 Nvidia GPU technology conference poster
http://on-demand.gputechconf.com/gtc/2014/poster/pdf/P4163_concurrent_VoIP_speech_decoding.pdf 

### How to build
In a Linux flavor, be sure to install *make*, *g++* and *nvidia-cuda-toolkit* package. Then clone or dowload code.
```
make
```

### How to run
Experiment may be performed through
```
./bin/test
```

### License
Code of g722 decoder is highy resembles to ITU reference implementation that should be trated as IEFT wishes.
All remaining part is MIT.

