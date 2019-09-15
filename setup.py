import os
import re
import sys
import glob
import shutil
import platform
import subprocess

from setuptools import setup, Extension, Command
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md")) as f:
    long_description = f.read()

BUILD_TEMP = os.path.expanduser("~/Desktop/vecpsf_build")


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # dealing with crappy Dropbox
        if ("(") in os.path.abspath(self.build_temp):
            self.build_temp = BUILD_TEMP
            extdir = os.path.join(self.build_temp, os.path.basename(extdir))

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
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

    user_options = []

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
                    raise ValueError("%s is not a path inside %s" % (_path, here))
                print("removing %s" % os.path.relpath(_path))
                if os.path.isdir(_path):
                    shutil.rmtree(_path)
                else:
                    os.remove(_path)

        if os.path.exists(BUILD_TEMP):
            print("removing %s" % BUILD_TEMP)
            shutil.rmtree(BUILD_TEMP)


class DeployCommand(Command):
    user_options = []
    repository_url = "https://test.pypi.org/legacy/"

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        self.run_command("test")
        self.run_command("bdist_wheel")
        subprocess.check_call(
            ["twine", "upload", "--repository-url", self.repository_url, "dist/*"]
        )
        self.run_command("clean")


setup(
    name="psfmodels",
    version="0.1.0",
    author="Talley Lambert",
    author_email="talley.lambert@gmail.com",
    license="GPL-3.0",
    url="https://github.com/tlambert03/PSFmodels-py",
    description="Scalar and vectorial models of the microscope point spread function (PSF).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("psfmodels")],
    cmdclass=dict(build_ext=CMakeBuild, clean=CleanCommand, deploy=DeployCommand),
    zip_safe=False,
    install_requires=["numpy"],
    packages=["psfmodels"],
)
