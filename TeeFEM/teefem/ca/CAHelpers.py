__doc__ = 'Code Aster helpers'

import sys
import os
import datetime
import getopt
import re
import time
import subprocess

global debug
debug = True
def dprint(*args, **kwds):
    global debug
    if debug:
        print args[0]

class File(object):
    fileid = 0
    def __init__(self):
        File.fileid += 1
        self.id = File.fileid

class Export(object):
    def __init__(self):
        self.actions = 'make_etude'
        self.debug = 'nodebug'
        self.lang = 'eng'
        self.mode = 'interactif'
        self.ncpus = 1 
        self.version = 'STA11.1' # testing
        self.memjeveux = 512.0
        self.mem_aster = 256.0
        self.tpmax = 5*3600 # => 5 hours
        self.tps_job = 30
        self.memjob = 10*1024*1024 # => 10 GB
        self.scriptdir = sys.path[0]
        self.scriptname = sys.argv[0]
        self.base = False
        self.basename = os.path.splitext(self.scriptname)[0]
        self.files = []
    @property
    def commfiles(self):
        return [file.name for file in self.files.values() if file.type == 'comm']
    def addfile(self,ffname=None,fname=None,ext='comm',type=None,opts='D',unite=1):
        f = File()
        if ffname:
            f.name = ffname
        elif fname:
            f.name = os.path.join(self.scriptdir, fname)
        else:
            f.name = os.path.join(self.scriptdir, '.'.join([self.basename,ext]))
        f.ext = ext
        f.type = type if type else ext
        f.unite = unite
        f.opts = opts
        self.files.append(f)
    def writetofile(self,fname=None):
        if fname:
            self.exportfile = os.path.join(self.scriptdir, fname)
        else:
            self.exportfile = os.path.join(self.scriptdir, '.'.join([self.basename,'export']))
        with open(self.exportfile, 'w') as fh:
            fh.write('# Export file created %s\n'%(datetime.datetime.now()))
            fh.write('P actions %s\n'%(self.actions))
            fh.write('P version %s\n'%(self.version))
            fh.write('P debug %s\n'%(self.debug))
            fh.write('P lang %s\n'%(self.lang))
            fh.write('P mode %s\n'%(self.mode))
            fh.write('P ncpus %s\n'%(int(self.ncpus)))
            fh.write('A memjeveux %s\n'%(int(self.memjeveux)))
            fh.write('A tpmax %s\n'%(int(self.tpmax)))
            fh.write('P memjob %s\n'%(int(self.memjob)))
            for file in self.files:
                p0 = 'R' if file.type == 'base' else 'F'
                fh.write('%s %s %s %s %s\n'%(p0, file.type, file.name, file.opts, str(file.unite)))
        return self.exportfile

class Study(object):
    def __init__(self, *args, **kwds):
#        print args
#        print kwds
        self.export = Export()
        self.name = kwds['name']
        self.description = kwds['description']
        if kwds.get('impr_base', False):
            basedir = kwds.get('basedir', False)
            if not basedir:
                self.export.addfile(fname='_base/%s.base'%(self.name), type='base', opts='R', unite=0)
            else:
                self.export.addfile(ffname='%s/%s.base'%(basedir,self.name), type='base', opts='R', unite=0)
        self.export.addfile(fname='%s.comm'%(self.name), type='comm', opts='D', unite=1)
        self.export.addfile(fname='%s.mess'%(self.name), type='mess', opts='R', unite=6)
        self.export.addfile(fname='results/%s.resu'%(self.name), type='resu', opts='R', unite=8)
        if kwds.has_key('usebase'):
            basedir = kwds.get('basedir', False)
            if not basedir:
                self.export.addfile(fname='_base/%s.base'%(kwds['usebase']), type='base', opts='D', unite=2)
            else:
                self.export.addfile(ffname='%s/%s.base'%(basedir,kwds['usebase']), type='base', opts='D', unite=2)
        if kwds.has_key('additionalfiles'):
            for f in kwds['additionalfiles']:
               self.export.addfile(**f)

    def run(self):
        self.exportfile = self.export.writetofile(fname='%s.export'%(self.name))
        retcode = subprocess.call(['/usr/bin/env', 'as_run', self.exportfile])
        return retcode

def caf(**kwds): return dict(kwds)

def usage(studies): 
    print("Studies")
    for s in studies:
        print("    {0:15s} {1}".format(s.name, s.description))

def run(studies, argv):
#    argv = sys.argv[1:]
#    print studies, argv
    study = {}
    for s in studies:
        study[s.name] = s.run
#    print argv
#    print study
    try:
        opts, args = getopt.getopt(argv, "dh", ['debug', 'help'] + study.keys())
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    for opt in opts:
        if opt in ("-d", "--debug"):
            global _debug
            _debug = True
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
    if len(args) == 0:
        usage(studies)
        sys.exit()
    for arg in args:
        if arg in study.keys():
            print "Running %s modelisation"%(arg) 
            study[arg]()

    
