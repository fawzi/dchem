---
title: Dchem
layout: Main
---
Dchem is a program package that

- is highly modular
- is written in a modern language, high level, yet easy to optimize down to the metal
- is parallel and can take advantage of etherogeneous clusters
- can easily be extended to use other programs to calculate energy and forces.

It allows you to explore the potential energy surface. It build a coarse incremental representation of the PES that is used to efficiently extend the exploration, and can be used for different kinds of analysis. The current version of the program is already quite general and can be applied to a large class of systems, in particular it can:

- identify minima and basins of attractions
- find which minima are directly connected through barriers that are smaller than the explored threshold
- subset a system to see only some degrees of freedom as active and to be explored
- add constraints to force a given alignement
- block the exploration of systems where particles "fly away"
- keep into account permutational symmetry.

This is the first release of the program, there might still be bugs.

Dchem needs the [blip](http://fawzi.github.com/blip) library, where it is also explained [how to setup](http://fawzi.github.com/blip/HowToD.html) a D environment, and getting started with the development, but to simply use it you can also try to use one of the  binaries that should appear here soon.

- [About](About.html) dchem
- [Starting out](Starting.html) with dchem
- The source code is available at [https://github.com/fawzi/dchem](https://github.com/fawzi/dchem)

News

- 21.2.2011- the dchem code goes public!
