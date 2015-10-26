#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 14:53:55 2015

@author: Sebastian Kelm (https://github.com/skelm)
@author: Rohith Mohan (https://github.com/rohithmohan)
@author: Willy Evangeslista (https://github.com/wevangelista)
"""

import sys
import os
import subprocess
import tempfile
import shutil
from glob import glob
from multiprocessing import Pool, cpu_count
from os.path import expanduser

HOME = expanduser("~")
CG2AT_path = HOME+"/MemProtMD/CG2AT/cg2at-nosolvate.pl"


def read_title(fname):
    with open(fname, "rb") as fin:
        for line in fin:
            if line.startswith("TITLE"):
                return line
    return ""


def read_pdb_change_title(fname, title=""):
    with open(fname, "rb") as fin:
        if not title:
            return fin.read()
        lines = []
        for line in fin:
            if line.startswith("TITLE"):
                lines.append(title)
            else:
                lines.append(line)
        return "".join(lines)

    
def write_file(fname, content):
    with open(fname, "wb") as fout:
        fout.write(str(content))


def run_cg2at(in_fname, out_fname):
    tempdir = tempfile.mkdtemp()
    try:
        fname = os.path.join(tempdir, "input.pdb")
        shutil.copy(in_fname, fname)
        p = subprocess.Popen(CG2AT_path+" input.pdb", cwd=tempdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        cg2at_outfile = ""
        for outf in ("atomistic-system.pdb", "atomistic-protein-lipid.pdb", "atomistic-protein.pdb", "atomistic-lipid.pdb"):
            outf = os.path.join(tempdir, "CG2AT", outf)
            if os.path.exists(outf):
                cg2at_outfile = outf
                break
        if cg2at_outfile:
            title = read_title(in_fname)
            with open(out_fname, "wb") as fout:
                #fout.write("MODEL\n")
                fout.write(read_pdb_change_title(cg2at_outfile, title))
                #fout.write("ENDMDL\n")
            #shutil.copy(cg2at_outfile, out_fname)
            shutil.rmtree(tempdir)
            return True
        else:
            print >>sys.stderr, "ERROR occurred while generating %s. Intermediate output files are here: %s" % (out_fname, tempdir)
            return False
    except:
    #except KeyboardInterrupt:
        shutil.rmtree(tempdir)
        raise
    #except:
    #    print >>sys.stderr, "ERROR occurred. Output files are here: "+tempdir
    #    #shutil.rmtree(tempdir)
    #    raise


def count_multi_pdb_models(infile):
    cryst1 = ""
    models = 0
    inmodel = False
    with open(infile) as f:
        for line in f:
            if line.startswith("CRYST1"):
                cryst1 = line
            elif line.startswith("ATOM  ") or line.startswith("HETATM"):
                inmodel = True
            elif line.startswith("ENDMDL"):
                if inmodel:
                    models += 1
                    inmodel = False
    if inmodel:
        models += 1
    assert cryst1, "Input file has no CRYST1 line. This is necessary for CG2AT to run."
    return models

def read_multi_pdb(infile, processes=1, processno=0):
    title = ""
    cryst1 = ""
    buffer = []
    modelno = 0
    skip = modelno % processes != processno
    with open(infile) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                if buffer:
                    yield (modelno, title+cryst1+"".join(buffer))
                    buffer = []
                modelno += 1
                skip = modelno % processes != processno
            elif not skip:
                if (line.startswith("ATOM  ") or line.startswith("HETATM")):
                    buffer.append(line)
                elif line.startswith("CRYST1"):
                    cryst1 = line
                elif line.startswith("TITLE"):
                    title = line
    if buffer:
        yield (modelno, title+cryst1+"".join(buffer))


def run_cg2at_batch(infile, outprefix, outsuffix_at="_at", outsuffix_cg="_cg", counter_format="%06d", delete_cg=False, verbose=False, resume=False, processes=1, processno=0):
    succeeded = 0
    total = 0
    for i, pdbtext in read_multi_pdb(infile, processes, processno):
        fname_in = outprefix + counter_format%i + outsuffix_cg + ".pdb"
        fname_out = outprefix + counter_format%i + outsuffix_at + ".pdb"
        if verbose:
            print "CG2AT", i, fname_out
        if resume and os.path.exists(fname_out):
            continue
        write_file(fname_in, pdbtext)
        success = run_cg2at(fname_in, fname_out)
        if success and delete_cg:
            os.remove(fname_in)
        if success:
            succeeded += 1
        total += 1
    return (succeeded, total)


def _run_cg2at_batch(args):
    try:
        return run_cg2at_batch(*args)
    except KeyboardInterrupt:
        raise RuntimeError("KeyboardInterrupt")


def run_cg2at_batch_parallel(infile, outprefix, outsuffix_at="_at", outsuffix_cg="_cg", counter_format="%06d", delete_cg=False, verbose=False, resume=False, processes=1):
    if processes == 1:
        return run_cg2at_batch(infile, outprefix, outsuffix_at, outsuffix_cg, counter_format, delete_cg, verbose, resume)
    
    if processes < 1:
        processes = cpu_count()
    n_models = count_multi_pdb_models(infile)
    if processes > n_models:
        processes = n_models
    
    if processes <= 1:
        return run_cg2at_batch(infile, outprefix, outsuffix_at, outsuffix_cg, counter_format, delete_cg, verbose, resume)
    
    pool = Pool(processes)
    args = []
    for processno in xrange(processes):
        args.append((infile, outprefix, outsuffix_at, outsuffix_cg, counter_format, delete_cg, verbose, resume, processes, processno))
    try:
        results = pool.map_async(_run_cg2at_batch, args).get()
    except KeyboardInterrupt:
        pool.terminate()
        raise
    
    succeeded = 0
    total = 0
    for s, t in results:
        succeeded += s
        total += t
    
    return (succeeded, total)
    

def concatenate_trajectory(infiles, outfile):
    p = subprocess.Popen("trjcat -f "+" ".join(infiles)+" -o "+outfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >>sys.stderr, "Failed to concatenate trajectory:\n"+out+err
    return p.returncode == 0


def convert_trajectory(in_trajectory, in_topology, out_trajectory, selection_center=1, selection_output=0, skip=0, pbc="", sep=False):
    "Converts a trajectory (e.g. XTC) to a different format (e.g. PDB)"
    STDIN = "%d\n%d\n" % (selection_center, selection_output)
    if pbc:
        pbc = "-pbc "+pbc
    if skip:
        skip = "-skip %d"%skip
    else:
        skip = ""
    if sep:
        sep = "-sep"
    else:
        sep = ""
    p = subprocess.Popen("gmx trjconv -f %s -s %s -o %s %s %s -center %s" % (in_trajectory, in_topology, out_trajectory, skip, pbc, sep), shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate(STDIN)
    
    if not os.path.exists(out_trajectory):
        print >>sys.stderr, "Failed to convert trajectory: %s"%(in_trajectory)
    
    
def remove_cg_solvent(pdbfile, outfile):
    "Remove the CG solvent and ion atoms from a PDB file"
    p = subprocess.Popen("grep -v ' W     W ' %s | grep -v ' ION ' > %s" % (pdbfile, outfile), shell=True)
    p.communicate()
    if not os.path.exists(outfile):
        print >>sys.stderr, "Failed to filter pdb file: %s"%(pdbfile)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='''This script takes a coarse-grained multi-model PDB file and converts
                       each model to all-atom representation using CG2AT. Default output is
                       a series of separate PDB files, but this can be changed via command
                       line options. Type --help for more info.''')
    parser.add_argument('--out_prefix', metavar="PREFIX", default="output", help='output file prefix (default: "output")')
    parser.add_argument('--out_suffix_at', metavar="SUFFIX", default="_at", help='output file suffix for atomistic representation (default: "_at")')
    parser.add_argument('--out_suffix_cg', metavar="SUFFIX", default="_cg", help='output file suffix for coarse-grained representation (default: "_cg")')
    parser.add_argument('--out_counter_format', metavar="FORMATSTRING", default="%06d", help='output file counter format (default: "%%06d")')
    #parser.add_argument('--delete_cg', default=False, action="store_true", help='conserve disk space by deleting coarse-grained intermediate output files once they are successfully converted to all-atom')
    parser.add_argument('--clean', default=False, action="store_true", help='conserve disk space by deleting all intermediate output files')
    parser.add_argument('--resume', default=False, action="store_true", help='enable resume mode: skip frames where output files already exist. This is meant to work when re-running the same command after aborting it with CTRL-C.')
    parser.add_argument('--verbose', default=False, action="store_true", help='print some status messages to STDOUT')
    parser.add_argument('--processes', metavar="CPUS", default=1, type=int, help='Number of processes to use. 0 or negative numbers mean using all CPUs. (default: 1)')
    parser.add_argument('--skip', metavar="N", default=0, type=int, help='output every Nth frame only')
    parser.add_argument('--center_protein', action="store_true", default=False, help='center on protein (default is to center the whole system)')
    parser.add_argument('--select_protein', action="store_true", default=False, help='select only the protein from the input trajectory (this is a lot faster)')
    parser.add_argument('trajectory', help='input trajectory file')
    parser.add_argument('topology', help='input topology file')
    args = parser.parse_args()
    
    delete_cg = True
    
    selection_center=0
    if args.center_protein:
        selection_center=1
    
    selection_output=0
    if args.select_protein:
        selection_output=1
    
    pdbtraj_temp = args.out_prefix+args.out_suffix_cg+"_solvated.pdb"
    pdbtraj = args.out_prefix+args.out_suffix_cg+".pdb"
    
    if not (args.resume and (os.path.exists(args.out_prefix+args.out_suffix_at) or os.path.exists(args.out_prefix+args.out_suffix_at+".xtc"))):
        
        if not (args.resume and os.path.exists(pdbtraj)):
            convert_trajectory(args.trajectory, args.topology, pdbtraj_temp, selection_center, selection_output, skip=args.skip, pbc="res")
            assert os.path.exists(pdbtraj_temp), "Failed to convert trajectory to multi-model PDB"
            
            remove_cg_solvent(pdbtraj_temp, pdbtraj)
            assert os.path.exists(pdbtraj), "Failed to filter out solvent from multi-model PDB file"
            
            os.remove(pdbtraj_temp)
        
        successes, total = run_cg2at_batch_parallel(pdbtraj, args.out_prefix, outsuffix_at=args.out_suffix_at, outsuffix_cg=args.out_suffix_cg, counter_format=args.out_counter_format, delete_cg=delete_cg, verbose=args.verbose, resume=args.resume, processes=args.processes)
        
        if args.verbose:
            print "%d of %d frames successfuly written" % (successes, total)
            print "concatenating trajectory using trjcat --> at.xtc and cg.xtc"
        
        shutil.copy(args.out_prefix+args.out_counter_format%0+args.out_suffix_at+".pdb", args.out_prefix+args.out_suffix_at+".pdb")
        
        concatenate_trajectory(glob(args.out_prefix+"*"+args.out_suffix_at+".pdb"), args.out_prefix+args.out_suffix_at)
        #concatenate_trajectory(glob(args.out_prefix+"*"+args.out_suffix_cg+".pdb"), args.out_prefix+args.out_suffix_cg)
    
    if os.path.exists(args.out_prefix+args.out_suffix_at) or os.path.exists(args.out_prefix+args.out_suffix_at+".xtc"):
        if args.clean:
            for fname in glob(args.out_prefix+"*"+args.out_suffix_at+".pdb"):
                if fname == args.out_prefix+args.out_suffix_at+".pdb":
                    continue
                os.remove(fname)
            if os.path.exists(pdbtraj):
                os.remove(pdbtraj)

