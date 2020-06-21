cmake_minimum_required(VERSION 3.12)
project(mycell)

set(CMAKE_CXX_STANDARD 17)

add_executable(mycell main.cpp Particle.cpp Particle.h ParticleType.cpp ParticleType.h System.cpp System.h Reaction.cpp Reaction.h json.h jsoncpp.cpp Interaction.cpp Interaction.h)