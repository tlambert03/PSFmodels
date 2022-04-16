import glob
import os
import platform
import re
import shutil
import subprocess
import sys
from distutils.version import LooseVersion

from setuptools import Command, Extension, setup
from setuptools.command.build_ext import build_ext

BUILD_TEMP = os.path.expanduser("~/Desktop/vecpsf_build")


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError as e:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            ) from e

        if platform.system() == "Windows":
            m = re.search(r"version\s*([\d.]+)", out.decode())
            assert m
            cmake_version = LooseVersion(m[1])
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # dealing with crappy Dropbox
        if ("(") in os.path.abspath(self.build_temp):  # type: ignore
            self.build_temp = BUILD_TEMP
            extdir = os.path.join(self.build_temp, os.path.basename(extdir))

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += [f"-DCMAKE_BUILD_TYPE={cfg}"]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        flg = env.get("CXXFLAGS", "")
        e = f'{flg} -DVERSION_INFO=\\"{self.distribution.get_version()}\\"'
        env["CXXFLAGS"] = e

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""

    CLEAN_FILES = ".eggs ./*.so ./build ./dist ./*.pyc ./*.tgz ./*.egg-info".split(" ")

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        here = os.path.abspath(os.path.curdir)
        for path_spec in self.CLEAN_FILES:
            # Make paths absolute and relative to this path
            abs_paths = glob.glob(os.path.normpath(os.path.join(here, path_spec)))
            for _path in [str(p) for p in abs_paths]:
                if not _path.startswith(here):
                    # Die if path in CLEAN_FILES is absolute + outside this directory
                    raise ValueError(f"{_path} is not a path inside {here}")
                print(f"removing {os.path.relpath(_path)}")
                if os.path.isdir(_path):
                    shutil.rmtree(_path)
                else:
                    os.remove(_path)

        if os.path.exists(BUILD_TEMP):
            print(f"removing {BUILD_TEMP}")
            shutil.rmtree(BUILD_TEMP)


setup(
    setup=[CMakeExtension("psfmodels")],
    cmdclass=dict(build_ext=CMakeBuild, clean=CleanCommand),
)
