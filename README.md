# Cammisol project

Cammisol is a Multi-Agent soil simulation project. Its purpose is to simulate the evolution of N/P/C cycles in soil, considering the heterogeneous nature of the soil and interactions between organisms and physical processes in the soil.

The Cammisol project is an experimentation workbench for expert and non-expert. Even if default parameters from our own case studies are provided, its purpose is to be empirically configured and calibrated to fit the wide diversity of possible environmental context. To do so, the Cammisol project is designed with a specific attention toward final users, ensuring that parameters of the model are meaningful and easy to set empirically to allow an easy testing of different hypothesis. To enhance the accessibility of the model to non-expert, Cammisol takes full advantage of the [Gama platform](https://gama-platform.org/) to provide explainable results in an user-friendly interface. Thanks to the high level and intuitive GAMA language (GAML), the source code of the model itself can be easily understood by people that are not computer scientists.

Considering the complexity of the simulated system, Cammisol should be seen as an high level framework that coordinates interactions between several expert submodels. For example, Cammisol provides an environment as a grid that contains C/N/P nutrients, where microbes colonies can live. The lifecycle and metabolism of each microbes colony is then managed from a dedicated submodel that contains C/N/P nutrients as input/output of the submodel. There is no need to know the details of the submodel from the high level Cammisol model. In consequence, several implementations of each submodel might be provided, in order to test different state of the art hypothesis.

# Installation

[GAMA 1.9.3](https://github.com/gama-platform/gama/releases/tag/1.9.3) is recommended to run the model.

The current project can then be imported directly in the Gama platform.

For more information about the Gama platform usage and [how to import the current project in Gama](https://gama-platform.org/wiki/ImportingModels), do not hesitate to check the [official Gama documentation](https://gama-platform.org/wiki/PlatformDocumentation).

# Structure of the project

The main project is divided in several subprojects:
- `cammisol`: contains all soil mechanics currently simulated by the model.
  - `environment`: features related to the grid structure and nutrients compartments of the soil.
  - `enzymatic_activity`: model dedicated to the decomposition of organic matter by enzymes, and to the optimisation of enzymatic activities according to objectives of a microbes colony.
  - `microbes`: features related to the metabolism of microbes colonies.
  - `nematode`: features related to the metabolism of nematodes.
  - `cammisol.gaml`: main Cammisol model, from which global experiments can be launched.
- `interface`: contains legacy code for a serious game interface. This project is currently **broken** and unmaintained: it should not be used until fixed. For reproducibility purposes, the [Cammisol 1.0 prerelease](https://github.com/u2worm/cammisol/releases/tag/v1.0-pre) might be used to run the game, using the associated Gama version. Notice however that the associated version of the model is obsolete, and do not include the most recent improvements.

The [Cammisol wiki](https://github.com/u2worm/cammisol/wiki) provides information about Cammisol models and useful experiments from a final user point of view. If more detailed explanations are required, do not hesitate to explore the code and the associated inline documentation using the Gama interface.

#Camisol model

The `cammisol` experiment of the `cammisol/cammisol.gaml` model allows to explore outputs of the model with default parameters.

Model parameters are dispatched in the global scope of model files.

# Participative interface

The participative version of the Camisol application can be launched from the `application.gaml` model, using the `application` experiment.

Notice that the participative version calls several instances of the standalone `cammisol` model (one per plot).
