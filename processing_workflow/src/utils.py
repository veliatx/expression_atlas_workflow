"""
Utility functions for processing workflow.
Contains functions for file integrity checks and other helper operations.
"""

from pathlib import Path 
import gzip 

def check_gzipped_fastq_integrity(str: fh) -> bool:
    """
    Check if a gzipped FASTQ file is valid and can be opened.
    
    Args:
        fh (str): File path to the gzipped FASTQ file
        
    Returns:
        (bool): True if file is valid and can be opened, False otherwise
    """
    if not fh:
        return True
    try:
        with gzip.open(fh, 'rb') as f:
            f.read(1) 
        return True 
    except (OSError, gzip.BadGzipFile) as e:
        print(f"Bad gzipped fastq: {fh},{e}")
        return False