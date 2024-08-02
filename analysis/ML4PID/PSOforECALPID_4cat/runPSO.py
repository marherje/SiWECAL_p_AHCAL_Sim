#!/usr/bin/env python
import argparse, os, subprocess

from PSO.PSOManager import PSOManager
from PSO.common     import *

if __name__ == '__main__':
    ### args
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', dest='config',
                        action='store', default='',
                        help='path to PSO configuration file')
    parser.add_argument('-o', '--output-dir', dest='output_dir',
                        action='store', default='',
                        help='path to output directory')
    parser.add_argument('--skip-trees', dest='skip_trees',
                        action='store_true', default=False,
                        help='skip creation of signal/background trees')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='print additional debug information')
    parser.add_argument('--dry-run', dest='dry_run',
                        action='store_true', default=False,
                        help='enable dry-run mode (PSO set up but not executed)')
    opts, opts_unknown = parser.parse_known_args()
    ###

    log_prx = os.path.basename(__file__)+' -- '

    if not opts.config    : KILL(log_prx+'unspecified path to PSO configuration file [-c]')
    if not opts.output_dir: KILL(log_prx+'unspecified path to output directory [-o]')

    if not os.path.isfile(opts.config): KILL(log_prx+'PSO configuration file not found [-c]: '+opts.config)

    if os.path.exists(opts.output_dir): KILL(log_prx+'path to output directory already exists [-o]: '+opts.output_dir)
    ###

    DATA_SUBDIR = 'InitData'

    print('Created dir '+opts.output_dir)

    OUTPUT_DIR = os.path.abspath(os.path.realpath(opts.output_dir))
    OUTPUT_SUBDIR = OUTPUT_DIR+'/'+DATA_SUBDIR

    subprocess.call(['mkdir', '-p', OUTPUT_DIR])
    subprocess.call(['mkdir', '-p', OUTPUT_SUBDIR])

    subprocess.call(['cp', opts.config, OUTPUT_SUBDIR+'/config.txt'])

    CONFIG_FPATH = OUTPUT_SUBDIR+'/config.txt'

    subprocess.call(['cp', '-r', 'PSO', OUTPUT_DIR])

    if not opts.skip_trees:
       subprocess.call(['root', '-b', '-q', OUTPUT_DIR+'/PSO/PrepareTrees.C+("'+CONFIG_FPATH+'", "'+OUTPUT_SUBDIR+'")'])

    PSO = PSOManager(opts.output_dir, DATA_SUBDIR, opts.verbose, CONFIG_FPATH)
    PSO.CompileAndSetupClientExecutable()
    PSO.InitParticles()

    if not opts.dry_run:
       PSO.RunPSO()
       PSO.PrintResult()
