# SubGradSolver

## Getting started

These instructions will let you build SubGradSolver.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp/wikis/custom).

Optionally, install the library in the system with:

```sh
sudo make install
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(SubGradSolver)
target_link_libraries(<my_target> SMS++::SubGradSolver)
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.
