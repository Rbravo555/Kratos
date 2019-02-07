<p align=center><img height="72.125%" width="72.125%" src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/kratos.png"></p>

[![Release][release-image]][releases] [![License][license-image]][license] [![Master][kratos-master-status]][travis-branches] [![appveyor-image]][appveyor-master]

_KRATOS Multiphysics_ ("Kratos") is a framework for building parallel, multi-disciplinary simulation software, aiming at modularity, extensibility, and high performance. Kratos is written in C++, and counts with an extensive Python interface. More in [Overview](https://github.com/KratosMultiphysics/Kratos/wiki/Overview)

**Kratos** is **free** under BSD-4 [license](https://github.com/KratosMultiphysics/Kratos/wiki/Licence) and can be used even in comercial softwares as it is. Many of its main applications are also free and BSD-4 licensed but each derived application can have its own propietary license.

[release-image]: https://img.shields.io/badge/release-6.0-green.svg?style=flat
[releases]: https://github.com/KratosMultiphysics/Kratos/releases

[license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/license.txt

[kratos-master-status]: https://travis-ci.org/KratosMultiphysics/Kratos.svg?branch=master
[travis-branches]: https://travis-ci.org/KratosMultiphysics/Kratos/branches

[appveyor-image]: https://ci.appveyor.com/api/projects/status/f9p57hci9ufkqkf5/branch/master?svg=true
[appveyor-master]: https://ci.appveyor.com/project/KratosMultiphysics/kratos


# Main Features
**Kratos** is __multiplatform__ and available for __Windows, Linux__ (several distros) and __macOS__.

**Kratos** is __OpenMP__ and __MPI__ parallel and scalable up to thousands of cores.

**Kratos** provides a core which defines the common framework and several application which work like plug-ins that can be extended in diverse fields.

# Contributions

In Kratos Core:
- [CIMNE](http://www.cimne.com) Fluid Dynamics
- [TUM](https://www.st.bgu.tum.de/) Lehrstuhl für Statik
- [Boost](http://www.boost.org/) for ublas
- [pybind11](https://github.com/pybind/pybind11) for exposing C++ to python
- [GidPost](https://www.gidhome.com/gid-plus/tools/476/gidpost/) providing output to [GiD](https://www.gidhome.com/)
- [AMGCL](https://github.com/ddemidov/amgcl) for its highly scalable multigrid solver
- [JSON](https://github.com/nlohmann/json) JSON for Modern C++
- [ZLib](https://zlib.net/) The compression library

In applications
- [Trilinos](https://trilinos.org/) for MPI linear algebra and solvers used in trilinos application
- [METIS](http://glaros.dtc.umn.edu/gkhome/views/metis) for partitioning in metis application
