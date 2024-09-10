from gzip import GzipFile
from typing import Any, BinaryIO, Generator, IO, TextIO, Union

FASTA_FILE_IGNORE_DELIMITERS = (b'>', b';')


class SequenceSegment:
    def __init__(self,
                 sequence_id: bytes,
                 data: bytes = b'',
                 epilogue: bool = False):
        self.id = sequence_id  # Sequence ID
        self.data = data  # Dinucleotide sequence byte string
        self.epilogue = epilogue  # Flag to mark end of sequence

    def is_empty(self):
        return len(self.data) == 0


def sequence_segments(
    fasta_file: Union[GzipFile, TextIO, IO[Any], BinaryIO],
    sequence_length: int,
    sequence_overlap_length: int = 0
) -> Generator[SequenceSegment, None, None]:

    """Iterates through a fasta file and yields SequenceSegment(s) for each
    sequence.

    The size of the each sequence segment is specified by the sequence_length
    parameter and will fill until there is sequence data in the fasta sequence
    left.
    The sequence lookahead length allows for subsequent iterations to have a
    length of segment that includes the previous sequence.
    """

    # NB: Immutable sequence of bytes
    current_sequence_id = b''
    # NB: Mutable sequence of bytes
    # NB: This is over a 1000x (not a typo) speed-up over a byte object
    working_sequence_buffer = bytearray()

    # NB: Line buffered reading is probably the best way to handle edge cases
    # around sequence ID parsing and newline characters
    # NB: fasta_file is assumed to be opened in binary mode (rb)
    for fasta_line in fasta_file:
        # If we are on a new sequence (sequence ID)
        # NB: Assume that either of the delimiters are indicators of a
        # new sequence, notably including comments
        if fasta_line.startswith(FASTA_FILE_IGNORE_DELIMITERS):  # type: ignore
            # Yield the current sequence segment if there is remaining sequence
            # NB: We always keep the lookahead/overlap in the working sequence
            # buffer, therefore there can only be sequence remaining if it is
            # longer than the lookahead/overlap length
            if len(working_sequence_buffer) > sequence_overlap_length:
                yield SequenceSegment(current_sequence_id,  # type: ignore
                                      bytes(working_sequence_buffer))

            # Get the new reference sequence name
            # NB: remove leading '>'
            current_sequence_id = fasta_line.split()[0][1:]  # type: ignore
            # Reset the working sequence buffer
            working_sequence_buffer = bytearray()

        # Otherwise the line we are on is sequence data
        else:
            fasta_line = fasta_line.rstrip()  # Remove trailing newline
            # Add to the working sequence buffer
            working_sequence_buffer += fasta_line  # type: ignore
            # While we have enough sequence buffer to fill a sequence segment
            while len(working_sequence_buffer) >= sequence_length:
                yield SequenceSegment(
                    current_sequence_id,  # type: ignore
                    bytes(working_sequence_buffer[:sequence_length]))
                # Truncate the working sequence buffer by the sequence length
                # minus the lookahead
                # XXX: Assert that the kmer/sequence length is always larger
                # than the lookahead length?
                truncate_length = sequence_length - sequence_overlap_length
                working_sequence_buffer = \
                    working_sequence_buffer[truncate_length:]

    # Yield the last sequence segment
    # NB: We always keep the lookahead in the working sequence buffer
    # So there needs to be check if it is longer the lookahead length
    if len(working_sequence_buffer) > sequence_overlap_length:
        yield SequenceSegment(current_sequence_id,  # type: ignore
                              bytes(working_sequence_buffer),
                              epilogue=True)
