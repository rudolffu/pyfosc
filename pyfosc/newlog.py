#!/usr/bin/env python
# coding: utf-8
import logging
import sys
from typing import Optional, Dict
from colorama import Fore, Back, Style


class ColoredFormatter(logging.Formatter):
    """Colored log formatter.
    https://gist.github.com/joshbode/58fac7ababc700f51e2a9ecdebe563ad
    """
    def __init__(self, *args, colors: Optional[Dict[str, str]]=None, **kwargs) -> None:
        """Initialize the formatter with specified format strings."""
        super().__init__(*args, **kwargs)
        self.colors = colors if colors else {}

    def format(self, record) -> str:
        """Format the specified record as text."""
        record.color = self.colors.get(record.levelname, '')
        record.reset = Style.RESET_ALL
        return super().format(record)

formatter = ColoredFormatter(
    '{color} [{levelname:5}] {reset} {asctime} | {name} | {message}',
    style='{', datefmt='%Y-%m-%dT%X%Z',
    colors={
        'DEBUG': Fore.CYAN,
        'INFO': Fore.GREEN,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.RED + Back.WHITE + Style.BRIGHT,
    }
)

def getLogger(name, level=logging.INFO):
    logger = logging.getLogger(name)
    logger.setLevel(level)    
    # Remove existing handlers
    if logger.hasHandlers():
        for handler in logger.handlers:
            logger.removeHandler(handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.set_name('console')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    return logger
