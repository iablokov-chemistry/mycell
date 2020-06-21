#!/bin/bash
g++-8 -Wall -o mycell -std=c++17 \
 ./main.cpp \
 ./Particle.h \
 ./Particle.cpp \
 ./ParticleType.h \
 ./ParticleType.cpp \
 ./System.h \
 ./System.cpp \
 ./Reaction.h \
 ./Reaction.cpp \
 ./Interaction.h \
 ./Interaction.cpp \
 ./json.h \
 ./jsoncpp.cpp
