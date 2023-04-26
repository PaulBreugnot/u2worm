# Camisol project

Camisol is a Multi-Agent soil simulation project.

# Directory structure

The main project is divided into two subprojects:
- `camisol`: contains all soil mechanics currently simulated by the model.
- `interface`: contains all UI features.

# Camisol model

The `camisol` experiment of the `camisol/camisol.gaml` model allows to explore
outputs of the model with default parameters.

Model parameters are dispatched in the global scope of model files.

# Participative interface

The participative version of the Camisol application can be launched from the
`application.gaml` model, using the `application` experiment.

Notice that the participative version calls several instances of the standalone
`camisol` model (one per plot).
