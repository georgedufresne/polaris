#!/usr/bin/env python3
"""Same sCMD from a code tree's default settings: argv code_dir setup out.npy"""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
from isochronetoolbox.scmd import Setup, context
S = Setup(**{k: v for k, v in json.load(open(sys.argv[2])).items() if k in Setup.__dataclass_fields__})
C = context(S)
cats = [C.model(10286, 14.1)[1](np.random.RandomState(S.crn_seed), 15.48, 0.08),
        C.model(10392, 13.3)[1](np.random.RandomState(S.crn_seed), 15.52, 0.065)]
np.save(sys.argv[3], np.vstack(cats)); print(sys.argv[1], np.vstack(cats).shape)
