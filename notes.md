# Notes

## Terminology

Definitions of terms as defined in Typogenetics.

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

Note that the AA duplet does not code for an amino acid. It is reserved as a "punctuation mark" to mean end of enzyme. Multiple amino acid sequences can be created from a single strand during translation.

## Folding

Each amino acid has the possibility of inducing a kink in the enzyme. r indicates a right turn in the enzyme, l indicates a left, and s indicates no kind induced and the enzyme remains straight at that amino acid.

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

## Analogy to Real Biology

DNA chains are made up of nucleotides. Each nucleotide is made up of 1. a phosphate group, 2. a ribose sugar, and 3. a base. In Typogenetics we drop the first two components of a nucleotide and our strands are just composed of bases.

Enzymes are one type of protein. All proteins in Typogenetics can actively operate on strands, so we just refer to all proteins in Typogenetics as enzymes.

Transcription is the process of turning DNA into mRNA, which is then read by ribosomes to create proteins through the process of translation. In Typogenetics we skip this step and our program (the player of the TNT game) does the job of a ribosome, creating enzymes without any machinery built from genetic code.

Proteins are actually made up of 20 distinct amino acids, compared to the 15 in Typogenetics.

Typically about 300 amino acids make up a complete protein. Strands of DNA can be made up of hundreds of thousands or millions of nucleotides. Compare to Typogenetics where strands are potentially much shorter and amino acid sequences are potentially on the order of the length of strands of bases.

Three consecutive bases/nucleotides form a "codon". Since there are 4 bases and each codon has 3 bases, the number of possible codons is 4x4x4 = 64. Seeing as there are 20 amino acids, multiple codons will map to a single amino acid.

Finally, in real biology there is no 1:1 relationship between an amino acid and some operation. The tertiary structure of a protein decides the function of the protein. It is the full context of the protein that determines how any one amino acid will function.

## Next Steps

Using a codon of more bases allows one strand to code for multiple proteins depending on the binding site. Increasing the codon size to 4 bases could increase the density of genes on a given strand. The instruction set can stay the same and contain redundancy as in actual biology.

Adding a third category of bases could potentially give the system more expressive power.
