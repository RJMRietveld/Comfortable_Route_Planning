This repository is part of the Bachelor End Project (BEP)

Comfortable  Route  Planning

from mechanical engineering of the Delft university of technology.

Authors: Bas de Jager, Maurits Neele, Mike Rietveld and Julian waas

## Structure of the code
your_route.py<br />
    -> route_generator_definitions.py<br />
    -> generate_least_sickening_route_definitions.py<br />
        --> route_generator_definitions.py<br />
        --> route_calculations_definitions.py<br />
        --> speed_and_acceleration_definitions.py<br />
        --> motion_sickness_definitions.py<br />
        --> LtiManip.py<br />

## Installations
Install [OSMNx](https://osmnx.readthedocs.io/en/stable/).<br />
You can install OSMnx with conda:<br />
```
conda config --prepend channels conda-forge
conda create -n ox --strict-channel-priority osmnx
```
