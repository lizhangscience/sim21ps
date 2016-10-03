#!/usr/bin/env python3
#
# Copyright (c) 2016 Zhixian MA <zxma_sjtu@qq.com>
# MIT license

"""
A test script to learn the grammar and code tyle of python, and try to make interesting docstrings.
"""

import os
import sys
import argparse
import logging

import numpy as np
from astropy.io import fits

import fg21sim
from fg21sim.configs import configs
from fg21sim.utils import setup_logging


