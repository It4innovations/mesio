
import sys, os, logging, subprocess, types

def configure(ctx):
    ctx.env.static = ctx.options.static
    ctx.link_cxx = types.MethodType(link_cxx, ctx)

    ctx.msg("Setting int width to", ctx.options.intwidth)
    ctx.msg("Setting build mode to", ctx.options.mode)

    """ Set compilers """
    ctx.find_program(ctx.options.mpicxx, var="MPICXX")
    ctx.load(ctx.options.cxx)
    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPICXX

    """ Set default compilers flags"""

    ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
    ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])

    if ctx.options.intwidth == "32":
        ctx.env.append_unique("DEFINES", [ "esint=int" ])
        ctx.env.append_unique("DEFINES_API", [ "MESIO_INT_WIDTH=32" ])
    if ctx.options.intwidth == "64":
        ctx.env.append_unique("DEFINES", [ "esint=long" ])
        ctx.env.append_unique("DEFINES_API", [ "MESIO_INT_WIDTH=64" ])

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++11", "-Wall" ])
    ctx.env.append_unique("CXXFLAGS", ctx.options.cxxflags.split())
    if ctx.options.mode == "release":
        ctx.env.append_unique("CXXFLAGS", [ "-O3", "-g" ])
    if ctx.options.mode == "devel":
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-g" ])
    if ctx.options.mode == "debug":
        ctx.env.append_unique("CXXFLAGS", [ "-O0", "-g" ])
    if ctx.options.mode == "profile":
        ctx.env.append_unique("CXXFLAGS", [ "-O3" ])

    ctx.env.append_unique("INCLUDES", "src")
    ctx.env.append_unique("DEFINES", [ "__ESMODE__="+ctx.options.mode.upper() ])

    """ Recurse to third party libraries wrappers"""
    recurse(ctx)
    print_available(ctx)

def build(ctx):
    def is_git_directory(path = '.'):
        return subprocess.call(['git', '-C', path, 'status'], stderr=subprocess.STDOUT, stdout = open(os.devnull, 'w')) == 0

    if is_git_directory():
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip().decode()
    else:
        commit = "unknown"

    ctx.env.append_unique("DEFINES_INFO", [ '__ESCOMMIT__=\"{0}\"'.format(commit) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXX__=\"{0}\"'.format(ctx.env.CXX[0]) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESBUILDPATH__=\"{0}\"'.format(ctx.bldnode.abspath()) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXXFLAGS__=\"{0}\"'.format(" ".join(ctx.env.CXXFLAGS)) ])

    # dirty hack
    # find better solution by waf
    ctx.env["STLIB_MARKER"] = ["-Wl,-Bstatic,--start-group"]
    ctx.env.prepend_value("SHLIB_MARKER", "-Wl,--end-group")

    features = "cxx cxxshlib"
    ctx.lib = ctx.shlib
    if ctx.env.static or ctx.options.static:
        features = "cxx"
        ctx.lib = ctx.stlib

    prefix = "mio"
    ctx.mesio = []

    def build(files, target, use=[]):
        ctx(features=features, source=files,target=prefix+target, use=use)
        ctx.mesio += [ prefix+target ] + use

    build(ctx.path.ant_glob('src/esinfo/**/*.cpp'), "esinfo", [ "INFO" ])
    build(ctx.path.ant_glob('src/config/**/*.cpp'), "config")
    build(ctx.path.ant_glob('src/basis/**/*.cpp'), "basis")
    build(ctx.path.ant_glob('src/wrappers/mpi/**/*.cpp'), "wmpi")

    build(ctx.path.ant_glob('src/mesh/**/*.cpp'), "mesh")
    build(ctx.path.ant_glob('src/input/**/*.cpp'), "input")
    build(ctx.path.ant_glob('src/output/**/*.cpp'), "output", [ "wpthread" ])
    build(ctx.path.ant_glob('src/wrappers/pthread/**/*.cpp'), "wpthread", [ "PTHREAD" ])
    build(ctx.path.ant_glob('src/wrappers/hdf5/**/*.cpp'), "whdf5", [ "HDF5" ])
    build(ctx.path.ant_glob('src/wrappers/metis/**/*.cpp'), "wmetis", [ "METIS" ])
    build(ctx.path.ant_glob('src/wrappers/parmetis/**/*.cpp'), "wparmetis", [ "PARMETIS" ])
    build(ctx.path.ant_glob('src/wrappers/scotch/**/*.cpp'), "wscotch", [ "SCOTCH" ])
    build(ctx.path.ant_glob('src/wrappers/ptscotch/**/*.cpp'), "wptscotch", [ "PTSCOTCH" ])
    build(ctx.path.ant_glob('src/wrappers/kahip/**/*.cpp'), "wkahip", [ "KAHIP" ])

    ctx.lib(source="src/api/wrapper.mesio.cpp", target="mesioapi", includes="include", use=ctx.mesio + ["API"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source=["src/api/api.mesio.cpp"], target="test.mesio", includes="include", use=ctx.mesio + ["API", "mesioapi"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source="src/app/mesio.cpp", target="mesio", use=ctx.mesio, stlib=ctx.options.stlibs, lib=ctx.options.libs)

def options(opt):
    opt.compiler = opt.add_option_group("Compiler options")
    opt.decomposers = opt.add_option_group("Third party graph partition tools")
    opt.other = opt.add_option_group("Other third party libraries")

    opt.compiler.add_option("--mpicxx",
        action="store",
        type="string",
        metavar="MPICXX",
        default=os.path.basename(os.getenv("MPICXX")) if os.getenv("MPICXX") else "mpic++",
        help="MPI compiler used for building of the library")

    opt.compiler.add_option("--cxx",
        action="store",
        type="string",
        metavar="CXX",
        default=os.path.basename(os.getenv("CXX")) if os.getenv("CXX") else "g++",
        help="C++ compiler")
    opt.compiler.add_option("--cxxflags",
        action="store",
        type="string",
        default="",
        help="C++ compiler flags (space separated list)")

    opt.compiler.add_option("--stlibs",
        action="store",
        type="string",
        default="",
        help="Additional static libraries")

    opt.compiler.add_option("--libs",
        action="store",
        type="string",
        default="",
        help="Additional dynamic libraries")

    opt.compiler.add_option("--intwidth",
        action="store",
        default="32",
        choices=["32", "64"],
        metavar="32,64",
        help="MESIO integer datatype width [default: %default]")

    modes=["release", "devel", "debug", "profile"]
    opt.compiler.add_option("-m", "--mode",
        action="store",
        default="release",
        choices=modes,
        help="MESIO build mode: " + ", ".join(modes) + " [default: %default]")

    opt.compiler.add_option("--static",
        action="store_true",
        default=False,
        help="MESIO executable file does not contain dynamic libraries.")

    recurse(opt)

def print_available(ctx):
    def _print(msg, err, libs, color="RED"):
        libs = [lib for lib in libs if "HAVE_" + lib.upper() in ctx.env["DEFINES_" + lib.upper()]]
        ctx.start_msg(msg)
        if len(libs) == 0:
            ctx.end_msg(err, color=color)
        else:
            ctx.end_msg("[ " + ", ".join(libs) + " ]", color="BLUE")
        return bool(len(libs))

    ctx.env["HAVE_DECOMPOSERS"] = _print(
        "Available graph partitioning tools",
        "not found [mesio cannot call a decomposer]",
        [ "metis", "parmetis", "scotch", "ptscotch", "kahip" ])
    _print(
        "Available miscellaneous libraries",
        "not found [mesio cannot load/store the XDMF format]",
        [ "pthread", "hdf5" ],
        "YELLOW")

""" Recurse to third party libraries wrappers"""
def recurse(ctx):
    """ MPI library """
    ctx.recurse("src/wrappers/mpi")

    """ Graph partition tools """
    ctx.recurse("src/wrappers/metis")
    ctx.recurse("src/wrappers/parmetis")
    ctx.recurse("src/wrappers/scotch")
    ctx.recurse("src/wrappers/ptscotch")
    ctx.recurse("src/wrappers/kahip")

    """ Other """
    ctx.recurse("src/wrappers/pthread")
    ctx.recurse("src/wrappers/hdf5")

from waflib import Logs
from waflib.Build import BuildContext
class ShowConfiguration(BuildContext):
    cmd = "show"
    fun = "show"

class ShowEnv(BuildContext):
    cmd = "env"
    fun = "env"

def show(ctx):
    ctx.logger = logging.getLogger('show')
    ctx.logger.handlers = Logs.log_handler()

    ctx.msg("CXX", " ".join(ctx.env.CXX))
    ctx.msg("DEFINES", " ".join(ctx.env.DEFINES))
    ctx.msg("CXXFLAGS", " ".join(ctx.env.CXXFLAGS))
    print_available(ctx)

def env(ctx):
    print(ctx.env)

def link_cxx(self, *k, **kw):
    includes = []
    libpath = []
    if "root" in kw and os.path.isdir(kw["root"]):
        includes = [ os.path.join(kw["root"], dir) for dir in os.listdir(kw["root"]) if dir.startswith("include") ]
        libpath = [ os.path.join(kw["root"], dir) for dir in os.listdir(kw["root"]) if dir.startswith("lib") ]

    general = dict(uselib_store=kw["name"].upper(), mandatory=False)
    if "mandatory" in kw:
        general["mandatory"] = kw["mandatory"]
    if "use" in kw:
        general["use"] = kw["use"]

    header = dict()
    self.env.stash()
    if "header_name" in kw:
        header = dict(header_name=kw["header_name"], define_name="", defines=["HAVE_" + kw["name"].upper()], includes=includes)
        header.update(general)
        header["msg"] = "Checking for '{0}' header".format(kw["name"])
        if not self.check_cxx(**header):
            self.env.revert()
            return False

        if "fragment" in kw:
            test = dict(execute=True)
            test.update(header)
            inc = [ "#include <{0}>\n".format(h) for h in kw["header_name"].split() ]
            test["fragment"] = "{0}int main(int argc, char** argv) {{ {1} }}".format("".join(inc), kw["fragment"])
            test["msg"] = "Checking for '{0}' settings".format(kw["name"])
            if not self.check_cxx(**test):
                self.env.revert()
                return False

    if "libs" in kw:
        libs = dict(stlib=kw["libs"], libpath=libpath, msg="Checking for '{0}' library".format(kw["name"]))
        libs.update(general)
        libs.update(header)
        if not self.options.static or not self.check_cxx(**libs):
            libs["lib"] = libs["stlib"]
            libs.pop("stlib")
            libs["msg"] = "Checking for '{0}' library".format(kw["name"])
            if not self.check_cxx(**libs):
                self.env.revert()
                return False

    return True
