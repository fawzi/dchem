---
title: Starting
layout: Default
---
Starting with Dchem
===================

To start using dchem there are two main thing that you need to understand:

- how to connect external programs
- and how the input of dchem is build

External programs
=================
To calculate energy and forces external programs can be used. Ideally one would like to be able to connect any program without modifying it, decouple the parallelization of the external program from the parallelization of the sampler, and minimize as much as possible the communication overhead between the two programs.
Dchem normally communicates with external programs either though files or through sockets.
Almost al programs have input files, can be started executing scripts and write output files.
Each calculator has a unique execution directory, and a shared template directory. When started the calculator copies all the files in the template directory to its execution directory. A file named `f.templ4`, is not copied, but generates the file `f` with the content of `f.templ4`, but where some tokens of the form `[keyword]` are replaced wherever the atoms change position.
The number at the end of the templ extension determines how often the file is updated: 

* 0: only the first time,
* 1: when the number or the kind of particles is changed,
* 2: when the position of the atom is changed,
* 3: when the position of the atoms is changed by a small amount,
* 4: when the position of the atoms is changed by a small amount in a smooth (extrapolable) way
* 5: when the configuration has not changed since last time

There are several special keywords, for example `[coordXyz]` is replaced by the coordinates in xyz format, `[coordTurbo]` by the coordinates in turbo mole format, `[evalE]` is 0 or one depending if the energy should be evaluated, `[changeLevel]` describes how much the positions have changed. The user can also define its own keywords.

The template approach allows one to create the inputs for his own program in a relatively simple way. To then calculate energies and forces it is possible to execute a command (that can for example call a script), and then parse some output files to collect energies and forces.
Turbomole and Molpro use this method and have a special parsing for the result.
Arbitrary programs can be connected if their output is prepared in one of the supported formats.

File based communication normally requires no modification of the external program, but is relatively slow. If the force evaluation is relatively quick, and a slow shared filesystem is used, file communication can become the bottleneck.
Dchem comes with a single c file that can be used to communicate arrays of basic types (reals, integers, chars) using sockets, form c or fortran programs.
Using it one can modify the external program implementing a simple communication protocol (set position, calculate energy/forces, get energy/forces) using sockets. This method can avoid the use of files, and can thus be much faster.
All the data being sent is serialized, something that is normally an acceptable overhead, and cleanly separates the parallelization of the external program.

The initial setup of cp2k uses files and templates, just like the other programs (needing basically no change to most of the program), but then the actual energy and force evaluation use the socket based communication, making updates to the configuration fast.
Such a strategy can be applied also to other programs to reduce the communication overhead with small changes to the external program.

Input
=====

dchem is a very modular program, and this is visible also in the input.
The input just describes various objects using a json-like syntax:

- strings without spaces can be written unquoted, any string can be written between double quotes ("bla \"bla")
- numbers are written in the normal format, fortran D specifier is accepted (1.3D-5)
- arrays are enclosed in square barckets ([1,2,3])
- dictionaries are enclosed in curly brackets and have the form
    { key1:value1, key2:value }
- objects just like dictionaries, but have class:className as first key
- The whole input is a dictionary, and objects in it can be referenced with
    { class:Ref, name:keyInTheInput }
- to simplify the writing most commas are optional and can be skipped

Objects can be anything, but there are some important kinds:

- samplers have a run method and do some operation
- methods devine an evaluation method, and a system on which it is applied
- configurations define a configuration
defining an object in the input creates it, and makes some minimal setup, but does not start its evaluation. In the input there is a special object called main, the object referenced by it is executed.

When doing a pnet based exploration (the normal way of exploring the potential energy surface with dchem) some other objects kinds are important:

- PNetSilos is a sampler, but is also the place where all the evaluatios are stored, it takes care of communicating topological infromation and its changes, and coordinate explorations. It refers several of the other object kinds
- silos workers perform an operation on the silos, they might load son data, or monitor and log some events. They are classified in normal workers and components that are registered with the silos and can be stopped. The silos has the loaders that are called before starting the exploration, the monitors that are called at regular intervals and the finishers that are called at the end of the exploration.
- work askers are samplers that connect to a silos and ask for work and evaluate energies and forces. A silos might start a work asker (evaluatorTask) automatically, or be executed alone (provided it has an url of a silos to connect to).
- explorers and observers are objects that register in the silos and get notified and can act during various important phases of the exploration. For example they might write a restart log (like PNetJournal). Explorers can also propose ne points to explore. Both these can be added to the explorers field of a silos.

Normally for an exploration your input will look like this:

- configuration
- method (that references the configuration)
- method limiter (to avoid having too many instances running concurrently)
- a loader to have some evaluations to start from
- loggers and explorers of pnet
- a work asker that actually computes things
- a silos that coordinates and starts everything
- main that references the silos

if you want to have just a worker then 

- configuration
- method (that references the configuration)
- method limiter (to avoid having too many instances running concurrently)
- a loader to have some evaluations to start from
- loggers and explorers of pnet
- a work asker that actually computes things
- a silos that coordinates and starts everything
- main that references the silos
