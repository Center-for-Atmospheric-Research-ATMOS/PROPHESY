TODO: uodate README

# ATTIRE

Framework for modeling kinetic energy analyzers

## Modular model

Should be implemented in the src folder some modules, one for each stage of an analyzer:

- lens
- hemispherical analyzer 
- detector: detection, electron multiplication, collector
- aperture (slit orbit mechanism)

The end user should probably not have to use the modules separatly (even though it should be possible) so, there should be a global function to call.

Also, it should contain one object per module that gather the main characteristic values to describe the modules.

## source, test and data description

TODO

- [test](./test/README.md)


## Refs

TODO: get all the bibliographic and library elements that can be used to understand what is done in the package, and also, some of the sources (e.g. Pyxel)
