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
    overlap_buffer = bytearray()

    sequences = []  # working list of sequences

    # import debugpy
    # debugpy.listen(5678)
    # debugpy.wait_for_client()
    # debugpy.breakpoint()
    # For every line in the fasta file
    for line in fasta_file:
        line = line.rstrip()  # Remove trailing newline

        # While there is enough working sequence buffer to fill the requested
        # sequence length
        # Create sequences for each segment
        sequences.extend(get_sequences_from_buffer(
                             working_sequence_buffer,
                             overlap_buffer,
                             sequence_length,
                             sequence_overlap_length))

        # If the current line is a sequence ID
        if line.startswith(FASTA_FILE_IGNORE_DELIMITERS):  # type: ignore
            # Yield the remaining sequences
            for sequence_buffer in get_remaining_sequence_segments(
                                   current_sequence_id,  # type: ignore
                                   sequences,
                                   working_sequence_buffer,
                                   overlap_buffer,
                                   sequence_length,
                                   sequence_overlap_length):
                yield sequence_buffer

            # Empty working sequences
            sequences.clear()

            # Update the working sequence ID
            # NB: remove leading '>'
            current_sequence_id = line.split()[0][1:]  # type: ignore
            # Clear the overlap buffer
            overlap_buffer.clear()
            # Clear the working buffer
            working_sequence_buffer.clear()

        # Otherwise the line is not a sequence ID and is sequence data
        else:
            # If any sequences were created
            if sequences:
                # Create all sequence segments but for the last
                for sequence in sequences[:-1]:
                    # Yield a sequence segment without the epilogue flag set
                    yield SequenceSegment(current_sequence_id, bytes(sequence))
                # Carry over the last sequence to the next iteration
                # in case this the last line it the sequence filled the
                # remaining buffer exactly
                sequences = [sequences[-1]]

            # Add the sequence line to the working buffer
            working_sequence_buffer += line  # type: ignore

    # Yield the remaining sequences
    for sequence_buffer in get_remaining_sequence_segments(
                           current_sequence_id,  # type: ignore
                           sequences,
                           working_sequence_buffer,
                           overlap_buffer,
                           sequence_length,
                           sequence_overlap_length):

        yield sequence_buffer


def get_sequences_from_buffer(working_sequence_buffer: bytearray,
                              overlap_sequence: bytearray,
                              sequence_length: int,
                              sequence_overlap_length: int):
    """Returns a list of overlapping byte sequences from a sequence buffer.
       Modifies the working sequence buffer and overlap buffer in place.
    """

    # If there is no sequence buffer
    if not working_sequence_buffer:
        # Return nothing
        return []

    sequences = []

    non_overlap_length = sequence_length - sequence_overlap_length

    while (len(working_sequence_buffer) + len(overlap_sequence) >=
           sequence_length):
        # If there is an overlap buffer
        if overlap_sequence:
            # Create the sequence with the overlap
            sequence = bytes(
                overlap_sequence +
                working_sequence_buffer[:non_overlap_length])
            bytes_used = non_overlap_length
        else:
            # Otherwise create the sequence without the overlap
            sequence = bytes(working_sequence_buffer[:sequence_length])
            bytes_used = sequence_length

        # Add to our working list of sequences
        sequences.append(sequence)
        # Update the overlap buffer if it exists by taking the last
        # current calculated sequence
        if sequence_overlap_length:
            # NB: Avoid re-assignment to modifiy in place
            overlap_sequence[:] = sequence[-sequence_overlap_length:]
        # Truncate the start of working sequence buffer by bytes used
        working_sequence_buffer[:bytes_used] = b''

    return sequences


def get_remaining_sequence_segments(sequence_id: bytes,
                                    sequences: list[bytes],
                                    working_sequence_buffer: bytearray,
                                    overlap_buffer: bytearray,
                                    sequence_length: int,
                                    sequence_overlap_length: int):
    # Assumes last of any buffer is the the epilogue

    sequence_segments = []

    sequences.extend(get_sequences_from_buffer(
                         working_sequence_buffer,
                         overlap_buffer,
                         sequence_length,
                         sequence_overlap_length))

    if working_sequence_buffer:
        sequences.append(bytes(overlap_buffer + working_sequence_buffer))

    # If any sequences were created
    if sequences:
        # Create a sequence segment for all but the last element
        # NB: Empty on a single list
        for sequence in sequences[:-1]:
            sequence_segments.append(
                    SequenceSegment(sequence_id, bytes(sequence))
            )

        # Create a sequence segment for the last element with the
        # epilogue flag set

        sequence_segments.append(
            SequenceSegment(sequence_id,
                            bytes(sequences[-1]),
                            epilogue=True)
        )

    return sequence_segments
