import logging
from enum import StrEnum, auto
from queue import Queue
from typing import Optional, Tuple

import numpy as np
from numpy.random import Generator

from typogenetics.typogenetics import Base, Rewriter, Strand, Translator

logger = logging.getLogger(__name__)


class EditType(StrEnum):
    MUTATE = auto()
    INSERT = auto()
    DELETE = auto()


class Editor:
    PROB_MUTATE = 0.80
    PROB_INSERT = 0.10
    PROB_DELETE = 0.10

    @classmethod
    def edit(cls, strand: Strand, rng: Generator) -> Strand:
        edit_type = cls.select_edit_type(rng)
        if edit_type == EditType.MUTATE:
            return cls.mutate(strand, rng)
        if edit_type == EditType.INSERT:
            return cls.insert(strand, rng)
        if edit_type == EditType.DELETE:
            return cls.delete(strand, rng)

    @classmethod
    def mutate(cls, strand: Strand, rng: Generator) -> Strand:
        r1 = rng.integers(0, len(strand))
        new_bases = strand.bases.copy()
        base = new_bases[r1]
        while new_bases[r1] == base:
            all_bases = [Base.A, Base.C, Base.G, Base.T]
            r2 = rng.integers(0, len(all_bases))
            new_bases[r1] = all_bases[r2]
        return Strand(new_bases)

    @classmethod
    def insert(cls, strand: Strand, rng: Generator) -> Strand:
        r1 = rng.integers(0, len(strand) + 1)
        new_bases = strand.bases.copy()
        all_bases = [Base.A, Base.C, Base.G, Base.T]
        r2 = rng.integers(0, len(all_bases))
        new_bases.insert(r1, all_bases[r2])
        return Strand(new_bases)

    @classmethod
    def delete(cls, strand: Strand, rng: Generator) -> Strand:
        r1 = rng.integers(0, len(strand))
        new_bases = strand.bases.copy()
        new_bases.pop(r1)
        return Strand(new_bases)

    @classmethod
    def select_edit_type(cls, rng: Generator) -> EditType:
        r = rng.random()
        edit_types = [
            (EditType.MUTATE, cls.PROB_MUTATE),
            (EditType.INSERT, cls.PROB_INSERT),
            (EditType.DELETE, cls.PROB_DELETE),
        ]
        assert np.isclose(sum(dict(edit_types).values()), 1.0)
        for edit_type, prob in edit_types:
            if r <= prob:
                return edit_type
            r -= prob
        raise ValueError("Random number is not in range [0, 1]")


class Search:
    @classmethod
    def random(
        cls,
        init_strand: Strand,
        n_iterations: int,
        rng: Generator,
        print_strands: bool = False,
    ) -> None:
        strands = [init_strand]
        known_set = set([str(init_strand)])
        for _ in range(n_iterations):
            enzyme_strand = strands[rng.integers(0, len(strands))]
            enzymes = Translator.translate(enzyme_strand)
            if len(enzymes) == 0:
                continue
            enzyme = enzymes[rng.integers(0, len(enzymes))]
            rewrite_strand = strands[rng.integers(0, len(strands))]
            new_strands = Rewriter.rewrite(enzyme, rewrite_strand)
            for strand in new_strands:
                if str(strand) not in known_set:
                    strands.append(strand)
                known_set.add(str(strand))

        if print_strands:
            print("Unique strands:")
            sorted_strands = sorted(known_set)
            for strand_str in sorted_strands:
                print(f"- {strand_str}")

        print(f"Discovered {len(known_set)} unique strands while simulating for {n_iterations} iterations")

    @classmethod
    def get_largest_rewrite_strand(
        cls,
        init_strand: Strand,
        apply_strand: Strand,
        log_rewrite: bool = False,
    ) -> Optional[Strand]:
        init_enzymes = Translator.translate(init_strand)
        if log_rewrite:
            logger.info(f"init_strand was translated into the enzymes: {init_enzymes}")
        init_enzyme_lengths = [len(enzyme) for enzyme in init_enzymes]
        if len(init_enzyme_lengths) == 0:
            return None
        largest_init_enzyme = init_enzymes[np.argmax(init_enzyme_lengths)]
        if log_rewrite:
            logger.info(f"The largest was {largest_init_enzyme}")

        strands = Rewriter.rewrite(largest_init_enzyme, apply_strand)
        if log_rewrite:
            logger.info(
                f"When the largest enzyme was applied to {apply_strand}, the following strands were produced: {strands}"
            )
        if len(strands) == 0:
            return None
        strand_lengths = [len(strand) for strand in strands]
        longest_strand = strands[np.argmax(strand_lengths)]
        if log_rewrite:
            logger.info(f"The longest strand produced was: {longest_strand}")
        return longest_strand

    @classmethod
    def bfs(
        cls,
        init_strand: Strand,
        apply_strand: Strand,
        target_depth: int,
        n_edits: int,
        rng: Generator,
        print_strands: bool = False,
    ) -> None:
        """Searches through the space of enzymes for enzymes with similar function.
        Every step of the search, edit a strand, translate it into enzymes, select the largest enzyme,
        apply it to the target strand to produce more strands, select the longest strand produced, compare that
        longest strand to the longest strand produced by the initial strand's corresponding largest enzyme.
        If the two strands match, then we assume that the enzyme has maintained its function after editing.
        """

        init_longest_rewrite_strand = cls.get_largest_rewrite_strand(init_strand, apply_strand, log_rewrite=True)
        if init_longest_rewrite_strand is None:
            logger.error("Could not find any rewrite strands for the given apply strand")
            return

        logger.info(
            f"Will search for enzymes that produce {init_longest_rewrite_strand} "
            f"when they are applied to {apply_strand}"
        )

        seen_strands = set()
        valid_strands = set()

        queue: Queue[Tuple[Strand, int]] = Queue()
        queue.put((init_strand, 0))
        while not queue.empty():
            curr_strand, depth = queue.get()

            if depth > target_depth:
                continue

            for _ in range(n_edits):
                edited_strand = Editor.edit(curr_strand, rng)
                if str(edited_strand) in seen_strands:
                    continue
                seen_strands.add(str(edited_strand))

                longest_rewrite_strand = cls.get_largest_rewrite_strand(edited_strand, apply_strand)
                if longest_rewrite_strand is None:
                    continue

                if longest_rewrite_strand == init_longest_rewrite_strand:
                    valid_strands.add(str(edited_strand))
                    queue.put((edited_strand, depth + 1))

        if print_strands:
            print("Valid strands:")
            sorted_strands = sorted(valid_strands)
            for strand_str in sorted_strands:
                print(f"- {strand_str}")

        print(
            f"Discovered {len(valid_strands)} valid strands while searching until "
            f"depth {target_depth} and branching factor {n_edits}"
        )
