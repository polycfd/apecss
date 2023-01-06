## How to contribute

### Reporting bugs
If you find a bug in APECSS, please report it on the [issue tracker](https://github.com/polycfd/apecss/issues/new?labels=bug).

### Suggesting enhancements
If you want to suggest a new feature or an improvement of a current feature, you can submit this
on the [issue tracker](https://github.com/polycfd/apecss/issues).

### Submitting a pull request
If you want to directly submit code to APECSS, you can do this by forking the repository, then submitting a pull request.

On opening a pull request, Github Actions conducts automated _build_ and _run_ tests. The _Build test_ checks whether the APECSS library compiles correctly on relevant operating systems. The _Run test_ is more comprehensive, in that it first compiles the APECSS library, then compiles and runs each [example](https://github.com/polycfd/apecss/tree/main/examples). Note that a successful _Run test_ does not imply correct results - the _Run test_ only tests the basic functionality of the code (e.g.~no segmentation faults), not the correctness of the results APECSS produces. 

### Code of conduct
We expect all our contributors to follow the [code of conduct](CODE_OF_CONDUCT.md). 

### Attribution
These contribution guidelines are adapted from the contribution guidelines of the [FEniCS/basix](https://github.com/FEniCS/basix) repository.
