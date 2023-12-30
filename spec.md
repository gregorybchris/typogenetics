# Typogenetics Specification

## Terminology

Definitions of terms as defined in Typogenetics:

| term                       | definition                                                              |
| -------------------------- | ----------------------------------------------------------------------- |
| bases                      | C, G, T, A                                                              |
| strand                     | string of bases                                                         |
| unit                       | position within a strand                                                |
| purines                    | A and G                                                                 |
| pyrimidines                | C and T                                                                 |
| complementary base pairing | when a strand is copied pyrimidines swap with purines, A <-> T, C <-> G |
| enzymes                    | operate on strands one unit at a time                                   |
| instruction                | an operation performed by an enzyme                                     |
| amino acid                 | three letter abbreviation for an instruction performed by an enzyme     |
| duplet                     | an adjacent pair of bases                                               |
| translation                | mapping from duplets to instructions                                    |
| primary structure          | amino acid sequence                                                     |
| tertiary structure         | folded structure of an enzyme                                           |
| gene                       | a portion of a strand that codes for a single enzyme                    |
| ribosome                   | reads strands and produces enzymes                                      |

## Instructions

Instructions as described in "Gödel, Escher, Bach":

| ins | action                                         |
| --- | ---------------------------------------------- |
| cut | cut strand(s)                                  |
| del | delete a base from strand                      |
| swi | switch enzyme to other strand                  |
| mvr | move one unit to the right                     |
| mvl | move one unit to the left                      |
| cop | turn on Copy mode                              |
| off | turn off Copy mode                             |
| ina | insert A to the right of this unit             |
| inc | insert C to the right of this unit             |
| ing | insert G to the right of this unit             |
| int | insert T to the right of this unit             |
| rpy | search for the nearest pyrimidine to the right |
| rpu | search for the nearest purine to the right     |
| lpy | search for the nearest pyrimidine to the left  |
| lpu | search for the nearest purine to the left      |

There are also some clarifications offered that preemptively address ambiguities that may have otherwise arisen:

- `cut` applies to both strands.
- `del` applies to only the strand on which the enzyme is working.
- `swi` moves the enzyme to the attached strand above the current strand. if there is no complementary base where the enzyme is currently bound, then the enzyme just detaches itself.
- insertion instructions will insert into both strands if Copy mode is on (with the complement inserted into the other strand). if Copy mode is off then a blank space is left in the complementary strand.
- if Copy mode is on and move or search instructions are encountered, then complementary bases should be manufactured everywhere the current strand slides.

## Translation

The first base from each duplet is on the y-axis and the second base is on the x-axis (ex. GC -> inc)

|     | A   | C   | G   | T   |
| --- | --- | --- | --- | --- |
| A   |     | cut | del | swi |
| C   | mvr | mvl | cop | off |
| G   | ina | inc | ing | int |
| T   | rpy | rpu | lpy | lpu |

> Note: The AA duplet does not code for an amino acid. It is reserved as a "punctuation mark" to mean end of enzyme. Multiple amino acid sequences can be created from a single strand during translation.

## Folding

Each amino acid has the possibility of inducing a kink in the enzyme. "r" indicates a right turn in the enzyme, "l" indicates a left turn, and "s" indicates no turn induced and the enzyme remains straight at that amino acid.

| ins | dir |
| --- | --- |
| cut | s   |
| del | s   |
| swi | r   |
| mvr | s   |
| mvl | s   |
| cop | r   |
| off | l   |
| ina | s   |
| inc | r   |
| ing | r   |
| int | l   |
| rpy | r   |
| rpu | l   |
| lpy | l   |
| lpu | l   |

## Binding-Preferences

The relative orientation of the first and last segments of an enzyme's tertiary structure determines the binding-preference of the enzyme.

Holding the orientation of the first segment to the right, the orientation of the last segment determines the binding-preference of the enzyme.

| first | last | base |
| ----- | ---- | ---- |
| R     | R    | A    |
| R     | U    | C    |
| R     | D    | G    |
| R     | L    | T    |

## Limitations to the Analogy

DNA chains are made up of nucleotides. Each nucleotide is made up of 1. a phosphate group, 2. a ribose sugar, and 3. a base. In Typogenetics we drop the first two components of a nucleotide and our strands are just composed of bases.

Enzymes are one type of protein. All proteins in Typogenetics can actively operate on strands, so we just refer to all proteins in Typogenetics as enzymes.

Transcription is the process of turning DNA into mRNA, which is then read by ribosomes to create proteins through the process of translation. In Typogenetics we skip this step and our program (the player of the Typogenetics game) does the job of a ribosome, creating enzymes without any machinery built from genetic code.

Proteins are actually made up of 20 distinct amino acids, compared to the 15 in Typogenetics.

Typically about 300 amino acids make up a complete protein. Strands of DNA can be made up of hundreds of thousands or millions of nucleotides. Compare to Typogenetics where strands are potentially much shorter and amino acid sequences are potentially on the order of the length of strands of bases.

Three consecutive bases/nucleotides form a "codon". Since there are 4 bases and each codon has 3 bases, the number of possible codons is 4x4x4 = 64. Seeing as there are 20 amino acids, multiple codons will map to a single amino acid.

Finally, in real biology there is no 1:1 relationship between an amino acid and some operation. The tertiary structure of a protein decides the function of the protein. It is the full context of the protein that determines how any one amino acid will function.

## Next Steps

### Larger codon size

If a codon is comprised of two nucleotides, then a single base can be translated into two different amino acids depending on the initial binding site. A codon of three nucleotides allows for three distinct meanings that a single base takes on, effectively increasing the density of genes without increasing the length of a strand. There may be a very good reason living systems on Earth use a codon of size three. I would be interested to explore the effects of codons of size 4, 5, 6. The size of the instruction set need not increase to accommodate an increased number of possible codons. Just as in real biology, a large diversity of nucleotide combinations can be mapped to a smaller set of amino acids with redundancy built in. Would increasing the density of genetic information on a strand help us evolve complex systems faster?

### More nucleobases

In real biology we have pyrimidines and purines. I would be curious to add a third category of bases. If I had to guess, C, T, G, A is close to the only code that satisfies both requirements of simplicity and error correction. Simplicity is a requirement because anything more complex would have been vanishingly unlikely to evolve out of primordial metabolic networks. And error correction, of course, to ensure genetic code would be stable enough to propagate over time. However we could have had a true binary code. Which makes me wonder what the effect would be of a hexadecimal code. Does increasing the number of available nucleotides increase the expressive power?

### Complex instruction set

The instruction set of 15 amino acids that Hofstadter gives us is certainly not the simplest possible instruction set, though there's something very beautiful about it being as reduced as it is. One does wonder how powerful strand rewriting could be with a few more instructions. I also wonder if the conditional rules that come in the box are a bit too complex, even. We currently can scan left/right until reaching a pyrimidine/purine. These are are conditionals, but not as simple as "if purine, move left one unit". Perhaps conditionals that simple could facilitate the evolution of more stable enzymes even if the enzymes need to be longer to do anything useful.

### Turing completeness

While I have not found anything definitive about whether Typogenetics is Turing complete, I would not be surprised if it were proven to be Turing incomplete. While there is certainly the ability to write to a tape, the lack of a set of states for the machine to be in is a bit worrying. Endowing an enzyme with a small finite state machine could be an interesting way to increase its representational power.

### Parallelism

Without changing the specification of Typogenetics at all it would be cool to speed up its execution by parallelizing. While each rewrite step is fundamentally serial, the processing of strands is an embarrassingly parallel operation. Especially if selection of strands and enzymes to interact is completely random, there are guaranteed to be no race conditions.

## Adding a spatial dimension

Inspired in part by Axelrod's experiments with agents playing the iterated prisoner's dilemma, you could limit strands to move around a "physical" space. Requiring interactions between enzymes and strands to be limited to spatially local interactions might promote more variation in evolved structures. More variation might come at the cost of lower complexity at first, but I can imagine some very improbable yet very destructive enzymes dominating if their radius of interaction is effectively infinite. Akin to ancient hydrothermal vents, rare pockets of fertile quiet may be necessary for fragile complexity to emerge slowly, undisturbed by its chaotic environment.
