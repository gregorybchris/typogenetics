# Extensions

## Larger codon size

If a codon is comprised of two nucleotides, then a single base can be translated into two different amino acids depending on the initial binding site. A codon of three nucleotides allows for three distinct meanings that a single base takes on, effectively increasing the density of genes without increasing the length of a strand. There may be a very good reason living systems on Earth use a codon of size three. I would be interested to explore the effects of codons of size 4, 5, 6. The size of the instruction set need not increase to accommodate an increased number of possible codons. Just as in real biology, a large diversity of nucleotide combinations can be mapped to a smaller set of amino acids with redundancy built in. Would increasing the density of genetic information on a strand help us evolve complex systems faster?

## More nucleobases

In real biology we have pyrimidines and purines. I would be curious to add a third category of bases. If I had to guess, C, T, G, A is close to the only code that satisfies both requirements of simplicity and error correction. Simplicity is a requirement because anything more complex would have been vanishingly unlikely to evolve out of primordial metabolic networks. And error correction, of course, to ensure genetic code would be stable enough to propagate over time. However we could have had a true binary code. Which makes me wonder what the effect would be of a hexadecimal code. Does increasing the number of available nucleotides increase the expressive power?

## Complex instruction set

The instruction set of 15 amino acids that Hofstadter gives us is certainly not the simplest possible instruction set, though there's something very beautiful about it being as reduced as it is. One does wonder how powerful strand rewriting could be with a few more instructions. I also wonder if the conditional rules that come in the box are a bit too complex, even. We currently can scan left/right until reaching a pyrimidine/purine. These are are conditionals, but not as simple as "if purine, move left one unit". Perhaps conditionals that simple could facilitate the evolution of more stable enzymes even if the enzymes need to be longer to do anything useful.

## Turing completeness

While I have not found anything definitive about whether Typogenetics is Turing complete, I would not be surprised if it were proven to be Turing incomplete. While there is certainly the ability to write to a tape, the lack of a set of states for the machine to be in is a bit worrying. Endowing an enzyme with a small finite state machine could be an interesting way to increase its representational power.

## Parallelism

Without changing the specification of Typogenetics at all it would be cool to speed up its execution by parallelizing. While each rewrite step is fundamentally serial, the processing of strands is an embarrassingly parallel operation. Especially if selection of strands and enzymes to interact is completely random, there are guaranteed to be no race conditions.

## Adding a spatial dimension

Inspired in part by Axelrod's experiments with agents playing the iterated prisoner's dilemma, you could limit strands to move around a "physical" space. Requiring interactions between enzymes and strands to be limited to spatially local interactions might promote more variation in evolved structures. More variation might come at the cost of lower complexity at first, but I can imagine some very improbable yet very destructive enzymes dominating if their radius of interaction is effectively infinite. Akin to ancient hydrothermal vents, rare pockets of fertile quiet may be necessary for fragile complexity to emerge slowly, undisturbed by its chaotic environment.

## Tuning

Many find it incredible that the John Horton Conway's Game of Life can produce and maintain so much complexity with such simple rules. I believe I remember Conway reacting to this impression in an interview once, saying something about how it's really not that incredible at all, given that the rules of the game were specifically selected in order to elicit that exact behavior of complexity and sustained complexity. I'm not sure if by that he meant that the rules were mathematically derived to produce the desired behavior or that the rules were tuned semi-blindly until the desired behavior emerged. Regardless, it has always intrigued me that if complex/interesting behavior does not initially emerge from a fairly complicated system, perhaps complex behaviors might emerge after fine-tuning parameters of that complicated system. Is there a way to parameterize the instructions of Typogenetics in such a way that they become tunable? Is there a metric we can optimize toward once we do have tunable instructions? If there's no good metric for complexity, what metric is worth optimizing?
