# Typogenetics

Implementation of the artificial evolutionary system described in Douglas Hofstadter's "Gödel, Escher, Bach".

Typogenetics (or "typographical genetics") is a toy version of real biological genetics. It simplifies a lot about how real genetics works to provide some useful intuitions with only a few moving parts.

**In computer terms**: You write a short program. The Typogenetics interpreter executes your program on some input you also provide and it produces some output. That output is guaranteed to be interpretable as another program. You can run that new program through the same interpreter. Given a good initial program, by running this execution in a loop we can hope to eventually produce a program that reproduces itself when interpreted.

**In biological terms**: Strands of bases are translated into sequences of amino acids. Amino acids are mapped to instructions. Instructions are applied to strands to produce new strands. This tangled hierarchy has the potential to enable the evolution of self-replicating enzymes.

For the nitty gritty of the Typogenetics specification see [spec.md](spec.md)

For some additional thoughts on ways to extend Typogenetics see [extensions.md](extensions.md)

## Requirements

- [Poetry](https://python-poetry.org/)

## Installation

```bash
poetry install
```

## Usage

```bash
# Translate a single strand into an enzyme
typo translate ATAGAGAGATCACATGTACGATAC

# Apply an enzyme to a strand to produce a set of new strands
typo rewrite cop-mvl-mvr-swi-cut-rpy AATACTAAACCGA

# Simulate many generations of evolution with a starting strand
typo simulate ATAGCGAATAGGATAATG --iter 10000 --seed 42

# Search for all strands that code for enzymes with similar function
typo search ATAAACGATAATTGACAGAGCGAATG ATCGATAGGGAACATGTCGT --edits 5 --depth 20 --seed 42
```

## Resources

### Repositories

- [andrejjocic/typogenetics](https://github.com/andrejjocic/typogenetics)
- [kortschak/typo](https://github.com/kortschak/typo)
- [lpalma/typogenetics](https://github.com/lpalma/typogenetics)
- [mkpankov/typogenetics](https://github.com/mkpankov/typogenetics)
- [nzni/typogenetics](https://github.com/nzni/typogenetics)
- [palm86/typogenetics](https://github.com/palm86/typogenetics)
- [Quuxplusone/TNT](https://github.com/Quuxplusone/TNT)

### Blogs

- [Gödel, Escher, Bach (David Fifield)](https://www.bamsoftware.com/hacks/geb)
- [A typogenetics implementation in Elixir (Danie Palm)](https://dev.to/palm86/a-typogenetics-implementation-in-elixir-1jfg)

### Papers

- [Autoreplicators and hypercycles in typogenetics (V. Kvasnicka, J. Pospichal)](https://www.sciencedirect.com/science/article/abs/pii/S016612800100464X)
- [Typogenetics: a logic of artificial propagating entities (Harold C. Morris)](https://open.library.ubc.ca/media/stream/pdf/831/1.0106810/1)
- [Typogenetics (Andrew Snare)](https://www.csse.monash.edu.au/hons/projects/1999/Andrew.Snare/thesis.pdf)
- [Typogenetics: an artificial genetic system (Louis Varetto)](https://pubmed.ncbi.nlm.nih.gov/8474250/)
