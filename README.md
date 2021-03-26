## About
A generic solver to solve moment models arise in kinetics gas theory.
Thanks to FEniCS, this solver is one-to-one in correspondence with 
mathematical formulation.
It is being developed as a part of PhD project **F2ME**. As the project
is currently undergoing code restructuring I have not made the 
documention of the code publically available (yet!). I plan to publish
it very soon. I am always interested in sharing information if you are trying
to use this solver. Drop me an email.

## Quick Usage
To make things simpler, I created a docker image which is in turn built
upon the FEniCS docker image. To get started, follow the instruction below:
Make sure you have `docker` and `docker-compose` installed.

```bash
$ git clone git@gitlab.com/19ec94/f2me.git
$ cd f2me
$ docker-compose run --rm fenics_rx_local
$ cd f2me/f2me 
$ python3 f2me.py input_example.yml
```

