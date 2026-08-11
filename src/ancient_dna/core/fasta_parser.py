"""
FASTA File Parser
=================

Reads standard FASTA format files and returns sequences as dictionaries.
Supports multi-sequence files and generator-based reading for large files.
"""

import logging
from typing import Dict, Generator, Optional, Tuple

logger = logging.getLogger(__name__)


def read_fasta_file(filename: str) -> Optional[Dict[str, str]]:
    """
    Read a FASTA file and return all sequences as a dictionary.

    Args:
        filename: Path to the FASTA file.

    Returns:
        Dictionary mapping sequence_name → sequence_data.
        Returns None if the file cannot be read.
    """
    sequences: Dict[str, str] = {}
    current_name: Optional[str] = None
    current_sequence: list = []

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Save previous sequence
                    if current_name:
                        sequences[current_name] = ''.join(current_sequence)
                    # Parse header — take only the name before first space for the key
                    current_name = line[1:].strip()
                    current_sequence = []
                elif current_name is not None:
                    current_sequence.append(line.upper())

            # Save the last sequence
            if current_name:
                sequences[current_name] = ''.join(current_sequence)

        logger.info("Parsed %d sequence(s) from '%s'", len(sequences), filename)
        return sequences

    except FileNotFoundError:
        logger.error("FASTA file not found: '%s'", filename)
        return None
    except Exception as e:
        logger.error("Error reading FASTA file '%s': %s", filename, e)
        return None


def read_fasta_generator(filename: str) -> Generator[Tuple[str, str], None, None]:
    """
    Read a FASTA file lazily using a generator — memory-efficient for large files.

    Yields:
        Tuples of (sequence_name, sequence_data).
    """
    current_name: Optional[str] = None
    current_sequence: list = []

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    if current_name:
                        yield current_name, ''.join(current_sequence)
                    current_name = line[1:].strip()
                    current_sequence = []
                elif current_name is not None:
                    current_sequence.append(line.upper())

            if current_name:
                yield current_name, ''.join(current_sequence)

    except FileNotFoundError:
        logger.error("FASTA file not found: '%s'", filename)
    except Exception as e:
        logger.error("Error reading FASTA file '%s': %s", filename, e)
