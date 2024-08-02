#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, subprocess

def KILL(log):
    raise SystemExit('\n '+'\033[1m'+'@@@ '+'\033[91m'+'FATAL'  +'\033[0m'+' -- '+log+'\n')
# --

def WARNING(log):
    print '\n '+'\033[1m'+'@@@ '+'\033[93m'+'WARNING'+'\033[0m'+' -- '+log+'\n'
# --

def EXE(cmd, suspend=True, verbose=False, dry_run=False):
    if verbose: print '\033[1m'+'>'+'\033[0m'+' '+cmd
    if dry_run: return

    _exitcode = os.system(cmd)

    if _exitcode and suspend: raise SystemExit(_exitcode)
# --

def get_output(cmd, permissive=False, warn=False):

    prc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if prc.returncode:

       log_msg = 'get_output -- shell command failed (execute command to reproduce the error):\n'+' '*14+'> '+cmd

       if not permissive:
          KILL(log_msg)

       elif warn:
          WARNING(log_msg)

       return None

    out, err = prc.communicate()

    return (out, err)
# --

def rreplace(str__, old__, new__, occurrence__):
    li_ = str__.rsplit(old__, occurrence__)
    return new__.join(li_)
# --

def which(program, permissive=False, warn=False):

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    exe_ls = []

    if fpath:
        if is_exe(program): exe_ls += [program]

    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)

            if is_exe(exe_file): exe_ls += [exe_file]

    if len(exe_ls) == 0:
        log_msg = 'which -- executable not found: '+program

        if permissive:
           if warn: WARNING(log_msg)
           return None

        else:
           KILL(log_msg)

    if len(exe_ls) >  1:
        log_msg = 'which -- executable "'+program+'" has multiple matches: \n'+str(exe_ls)

        if permissive:
           if warn: WARNING(log_msg)
           return None

        else:
           KILL(log_msg)

    return exe_ls[0]
# --

def is_int(value):

    try: int(value)
    except ValueError: return False

    return True
# --

def is_float(value):

    try: float(value)
    except ValueError: return False

    return True
# --

def HTCondor_jobIDs(username=None, permissive=False, warn=False):

    if not username:
       if 'USER' in os.environ: username = os.environ['USER']

    if not username:
       KILL('HTCondor_jobIDs -- undefined argument "username"')

    _condorq_jobIDs = {}

    _condorq_ret = get_output('condor_q '+str(username)+' -nobatch', permissive, warn)

    if _condorq_ret == None: return None

    _condorq_lines = _condorq_ret[0].split('\n')

    for _i_condorq in _condorq_lines:

        _i_condorq_pieces = _i_condorq.split()

        if len(_i_condorq_pieces) == 9:

           if is_float(_i_condorq_pieces[0]) and str(_i_condorq_pieces[1]) == username:

              _condorq_jobIDs[_i_condorq_pieces[0]] = {
                'ID'        : _i_condorq_pieces[0],
                'OWNER'     : _i_condorq_pieces[1],
                'SUBMITTED' : _i_condorq_pieces[2]+' '+_i_condorq_pieces[3],
                'RUN_TIME'  : _i_condorq_pieces[4],
                'STATUS'    : _i_condorq_pieces[5],
                'PRIORITY'  : _i_condorq_pieces[6],
                'SIZE'      : _i_condorq_pieces[7],
                'CMD'       : _i_condorq_pieces[8],
              }

    return _condorq_jobIDs

def HTCondor_job_executables(username=None, permissive=False, warn=False):

    if not username:
       if 'USER' in os.environ: username = os.environ['USER']

    if not username:
       KILL('HTCondor_jobIDs -- undefined argument "username"')

    _condorq_ret = get_output('condor_q '+str(username)+' -nobatch -long | grep "Cmd = "', permissive, warn)

    if _condorq_ret == None: return None

    _condorq_cmds = _condorq_ret[0].split('\n')

    _ret_paths = []

    for _i_cmd in _condorq_cmds:

        if _i_cmd == '': continue

        _i_cmd_pieces = _i_cmd.split(' = ')
        if len(_i_cmd_pieces) != 2: return None

        _exe_path = _i_cmd_pieces[1]
        _exe_path = _exe_path.replace(' ', '')

        if _exe_path.startswith('"'): _exe_path = _exe_path[+1:]
        if _exe_path.endswith  ('"'): _exe_path = _exe_path[:-1]

        _ret_paths += [os.path.abspath(os.path.realpath(_exe_path))]

    return _ret_paths

def HTCondor_executable_from_jobID(jobID):

    _condorq_ret = get_output('condor_q '+str(jobID)+' -long | grep "Cmd = "', permissive=True)

    if _condorq_ret == None: return None

    _condorq_cmd = _condorq_ret[0].split('\n')

    _condorq_cmd = [_tmp for _tmp in _condorq_cmd if _tmp != '']

    if len(_condorq_cmd) != 1: return None

    _condorq_cmd_pieces = _condorq_cmd[0].split(' = ')
    if len(_condorq_cmd_pieces) != 2: return None

    _exe_path = _condorq_cmd_pieces[1]
    _exe_path = _exe_path.replace(' ', '')

    if _exe_path.startswith('"'): _exe_path = _exe_path[+1:]
    if _exe_path.endswith  ('"'): _exe_path = _exe_path[:-1]

    _exe_path = os.path.abspath(os.path.realpath(_exe_path))

    return _exe_path
