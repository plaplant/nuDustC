from spack.package import *

class Nudustc(CMakePackage):
    """A model for dust production and evolution in astrophysical hydrocodes"""

    homepage = "https://github.com/lanl/nudustc/main/index.html"
    git = "https://github.com/lanl/nudustc.git"

    maintainers = ["mauneyc-LANL"]

    version("main", branch="main")

    variant("openmp", default=False, description="OpenMP support")
    variant("mpi", default=False, description="MPI support")

    depends_on("cmake@3.20:")
    depends_on("sundials@6.6.1: ~ARKODE~CVODES~IDA~IDAS~KINSOL cxxstd=17 ~examples~examples-install+int64")
    depends_on("boost@1.80: +program_options+filesystem+serialization cxxstd=17")

    depends_on("sundials+mpi", when="+mpi")
    depends_on("sundials~mpi", when="~mpi")
    depends_on("mpi", when="+mpi")

    def cmake_args(self):
        args = [
            self.define_from_variant("NUDUSTC_ENABLE_OPENMP", "openmp"),
            self.define_from_variant("NUDUSTC_ENABLE_MPI", "mpi")
        ]
        return args
    
