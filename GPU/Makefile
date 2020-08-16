CC = nvcc
SRC_EXT = .cu
OUT_EXT = .out
CFLAGS  = -O2
LIBES   = $(shell pkg-config --cflags --libs magma | sed -e 's/-fopenmp/-Xcompiler "-fopenmp"/g')
HEADER_DIR = ../Headers

ALL_O  = $(wildcard ${HEADER_DIR}/*.o)
ALL_H  = $(wildcard ${HEADER_DIR}/*.h)

TARGET  = testing_cheevd_GPU.out testing_zheevd_GPU.out
EXECS   = ${TARGET:.out=${OUT_EXT}}
SOURCES = ${TARGET:.out=${SRC_EXT}}

.SUFFIXES: .cu

all : ${FRC}
	${MAKE} "CFLAGS=${CFLAGS}" "FRC=${FRC}" ${EXECS}

Makefile : ${FRC}
	cp basic.mk $@
	@echo '# Automatically-generated dependencies list:' >> $@
	@${CC} -MM ${SOURCES} ${CFLAGS} ${LIBES}>> $@

.cu.out:
	${CC} -o $@ $< ${CFLAGS} ${ALL_O} ${LIBES}

force_rebuild :
# Automatically-generated dependencies list:
testing_cheevd_GPU.o : testing_cheevd_GPU.cu \
    /usr/local/cuda/include/cuda_runtime.h \
    /usr/local/cuda/include/crt/host_config.h \
    /usr/local/cuda/include/builtin_types.h \
    /usr/local/cuda/include/device_types.h \
    /usr/local/cuda/include/crt/host_defines.h \
    /usr/local/cuda/include/driver_types.h \
    /usr/local/cuda/include/vector_types.h \
    /usr/local/cuda/include/surface_types.h \
    /usr/local/cuda/include/texture_types.h \
    /usr/local/cuda/include/library_types.h \
    /usr/local/cuda/include/channel_descriptor.h \
    /usr/local/cuda/include/cuda_runtime_api.h \
    /usr/local/cuda/include/cuda_device_runtime_api.h \
    /usr/local/cuda/include/driver_functions.h \
    /usr/local/cuda/include/vector_functions.h \
    /usr/local/cuda/include/vector_functions.hpp \
    /usr/local/cuda/include/crt/common_functions.h \
    /usr/include/c++/7/new \
    /usr/include/x86_64-linux-gnu/c++/7/bits/c++config.h \
    /usr/include/x86_64-linux-gnu/c++/7/bits/os_defines.h \
    /usr/include/x86_64-linux-gnu/c++/7/bits/cpu_defines.h \
    /usr/include/c++/7/exception \
    /usr/include/c++/7/bits/exception.h \
    /usr/include/c++/7/bits/exception_ptr.h \
    /usr/include/c++/7/bits/exception_defines.h \
    /usr/include/c++/7/bits/cxxabi_init_exception.h \
    /usr/include/c++/7/typeinfo \
    /usr/include/c++/7/bits/hash_bytes.h \
    /usr/include/c++/7/bits/nested_exception.h \
    /usr/include/c++/7/bits/move.h \
    /usr/include/c++/7/bits/concept_check.h \
    /usr/include/c++/7/type_traits \
    /usr/include/c++/7/stdlib.h \
    /usr/include/c++/7/cstdlib \
    /usr/include/c++/7/bits/std_abs.h \
    /usr/local/cuda/include/crt/math_functions.h \
    /usr/include/c++/7/math.h \
    /usr/include/c++/7/cmath \
    /usr/include/c++/7/bits/cpp_type_traits.h \
    /usr/include/c++/7/ext/type_traits.h \
    /usr/local/cuda/include/crt/math_functions.hpp \
    /usr/local/cuda/include/cuda_surface_types.h \
    /usr/local/cuda/include/cuda_texture_types.h \
    /usr/local/cuda/include/crt/device_functions.h \
    /usr/local/cuda/include/crt/device_functions.hpp \
    /usr/local/cuda/include/device_atomic_functions.h \
    /usr/local/cuda/include/device_atomic_functions.hpp \
    /usr/local/cuda/include/crt/device_double_functions.h \
    /usr/local/cuda/include/crt/device_double_functions.hpp \
    /usr/local/cuda/include/sm_20_atomic_functions.h \
    /usr/local/cuda/include/sm_20_atomic_functions.hpp \
    /usr/local/cuda/include/sm_32_atomic_functions.h \
    /usr/local/cuda/include/sm_32_atomic_functions.hpp \
    /usr/local/cuda/include/sm_35_atomic_functions.h \
    /usr/local/cuda/include/sm_60_atomic_functions.h \
    /usr/local/cuda/include/sm_60_atomic_functions.hpp \
    /usr/local/cuda/include/sm_20_intrinsics.h \
    /usr/local/cuda/include/sm_20_intrinsics.hpp \
    /usr/local/cuda/include/sm_30_intrinsics.h \
    /usr/local/cuda/include/sm_30_intrinsics.hpp \
    /usr/local/cuda/include/sm_32_intrinsics.h \
    /usr/local/cuda/include/sm_32_intrinsics.hpp \
    /usr/local/cuda/include/sm_35_intrinsics.h \
    /usr/local/cuda/include/sm_61_intrinsics.h \
    /usr/local/cuda/include/sm_61_intrinsics.hpp \
    /usr/local/cuda/include/crt/sm_70_rt.h \
    /usr/local/cuda/include/crt/sm_70_rt.hpp \
    /usr/local/cuda/include/surface_functions.h \
    /usr/local/cuda/include/texture_fetch_functions.h \
    /usr/local/cuda/include/texture_indirect_functions.h \
    /usr/local/cuda/include/surface_indirect_functions.h \
    /usr/local/cuda/include/device_launch_parameters.h \
    /home/sugimoto/magma/include/magma_v2.h \
    /home/sugimoto/magma/include/magmablas.h \
    /home/sugimoto/magma/include/magma_copy.h \
    /home/sugimoto/magma/include/magma_types.h \
    /usr/local/cuda/include/cuda.h \
    /usr/include/c++/7/bits/stl_relops.h \
    /usr/include/c++/7/bits/stl_pair.h \
    /usr/include/c++/7/initializer_list \
    /usr/local/cuda/include/cusparse_v2.h \
    /usr/local/cuda/include/cusparse.h \
    /home/sugimoto/magma/include/magmablas_z.h \
    /home/sugimoto/magma/include/magmablas_c.h \
    /home/sugimoto/magma/include/magmablas_d.h \
    /home/sugimoto/magma/include/magmablas_s.h \
    /home/sugimoto/magma/include/magmablas_zc.h \
    /home/sugimoto/magma/include/magmablas_ds.h \
    /home/sugimoto/magma/include/magmablas_h.h \
    /home/sugimoto/magma/include/magma_batched.h \
    /home/sugimoto/magma/include/magma_zbatched.h \
    /home/sugimoto/magma/include/magma_cbatched.h \
    /home/sugimoto/magma/include/magma_dbatched.h \
    /home/sugimoto/magma/include/magma_sbatched.h \
    /home/sugimoto/magma/include/magma_hbatched.h \
    /home/sugimoto/magma/include/magma_vbatched.h \
    /home/sugimoto/magma/include/magma_zvbatched.h \
    /home/sugimoto/magma/include/magma_cvbatched.h \
    /home/sugimoto/magma/include/magma_dvbatched.h \
    /home/sugimoto/magma/include/magma_svbatched.h \
    /home/sugimoto/magma/include/magma_bulge.h \
    /home/sugimoto/magma/include/magma_zbulge.h \
    /home/sugimoto/magma/include/magma_cbulge.h \
    /home/sugimoto/magma/include/magma_dbulge.h \
    /home/sugimoto/magma/include/magma_sbulge.h \
    /home/sugimoto/magma/include/magma_z.h \
    /home/sugimoto/magma/include/magma_zgehrd_m.h \
    /home/sugimoto/magma/include/magma_c.h \
    /home/sugimoto/magma/include/magma_cgehrd_m.h \
    /home/sugimoto/magma/include/magma_d.h \
    /home/sugimoto/magma/include/magma_dgehrd_m.h \
    /home/sugimoto/magma/include/magma_s.h \
    /home/sugimoto/magma/include/magma_sgehrd_m.h \
    /home/sugimoto/magma/include/magma_zc.h \
    /home/sugimoto/magma/include/magma_ds.h \
    /home/sugimoto/magma/include/magma_auxiliary.h \
    /home/sugimoto/magma/include/magma_htc.h \
    /home/sugimoto/magma/include/magma_lapack.h \
    /home/sugimoto/magma/include/magma_mangling.h \
    /home/sugimoto/magma/include/magma_mangling_cmake.h \
    /home/sugimoto/magma/include/magma_zlapack.h \
    /home/sugimoto/magma/include/magma_clapack.h \
    /home/sugimoto/magma/include/magma_dlapack.h \
    /home/sugimoto/magma/include/magma_slapack.h \
    /usr/local/cuda/include/curand_kernel.h \
    /usr/local/cuda/include/curand.h \
    /usr/local/cuda/include/curand_discrete.h \
    /usr/local/cuda/include/curand_precalc.h \
    /usr/local/cuda/include/curand_mrg32k3a.h \
    /usr/local/cuda/include/curand_mtgp32_kernel.h \
    /usr/local/cuda/include/curand_mtgp32.h \
    /usr/local/cuda/include/curand_philox4x32_x.h \
    /usr/local/cuda/include/curand_globals.h \
    /usr/local/cuda/include/curand_uniform.h \
    /usr/local/cuda/include/curand_normal.h \
    /usr/local/cuda/include/curand_normal_static.h \
    /usr/local/cuda/include/curand_lognormal.h \
    /usr/local/cuda/include/curand_poisson.h \
    /usr/local/cuda/include/curand_discrete2.h \
    /usr/local/cuda/include/curand_mtgp32_host.h \
    /usr/local/cuda/include/curand_mtgp32dc_p_11213.h
testing_zheevd_GPU.o : testing_zheevd_GPU.cu \
    /usr/local/cuda/include/cuda_runtime.h \
    /usr/local/cuda/include/crt/host_config.h \
    /usr/local/cuda/include/builtin_types.h \
    /usr/local/cuda/include/device_types.h \
    /usr/local/cuda/include/crt/host_defines.h \
    /usr/local/cuda/include/driver_types.h \
    /usr/local/cuda/include/vector_types.h \
    /usr/local/cuda/include/surface_types.h \
    /usr/local/cuda/include/texture_types.h \
    /usr/local/cuda/include/library_types.h \
    /usr/local/cuda/include/channel_descriptor.h \
    /usr/local/cuda/include/cuda_runtime_api.h \
    /usr/local/cuda/include/cuda_device_runtime_api.h \
    /usr/local/cuda/include/driver_functions.h \
    /usr/local/cuda/include/vector_functions.h \
    /usr/local/cuda/include/vector_functions.hpp \
    /usr/local/cuda/include/crt/common_functions.h \
    /usr/include/c++/7/new \
    /usr/include/x86_64-linux-gnu/c++/7/bits/c++config.h \
    /usr/include/x86_64-linux-gnu/c++/7/bits/os_defines.h \
    /usr/include/x86_64-linux-gnu/c++/7/bits/cpu_defines.h \
    /usr/include/c++/7/exception \
    /usr/include/c++/7/bits/exception.h \
    /usr/include/c++/7/bits/exception_ptr.h \
    /usr/include/c++/7/bits/exception_defines.h \
    /usr/include/c++/7/bits/cxxabi_init_exception.h \
    /usr/include/c++/7/typeinfo \
    /usr/include/c++/7/bits/hash_bytes.h \
    /usr/include/c++/7/bits/nested_exception.h \
    /usr/include/c++/7/bits/move.h \
    /usr/include/c++/7/bits/concept_check.h \
    /usr/include/c++/7/type_traits \
    /usr/include/c++/7/stdlib.h \
    /usr/include/c++/7/cstdlib \
    /usr/include/c++/7/bits/std_abs.h \
    /usr/local/cuda/include/crt/math_functions.h \
    /usr/include/c++/7/math.h \
    /usr/include/c++/7/cmath \
    /usr/include/c++/7/bits/cpp_type_traits.h \
    /usr/include/c++/7/ext/type_traits.h \
    /usr/local/cuda/include/crt/math_functions.hpp \
    /usr/local/cuda/include/cuda_surface_types.h \
    /usr/local/cuda/include/cuda_texture_types.h \
    /usr/local/cuda/include/crt/device_functions.h \
    /usr/local/cuda/include/crt/device_functions.hpp \
    /usr/local/cuda/include/device_atomic_functions.h \
    /usr/local/cuda/include/device_atomic_functions.hpp \
    /usr/local/cuda/include/crt/device_double_functions.h \
    /usr/local/cuda/include/crt/device_double_functions.hpp \
    /usr/local/cuda/include/sm_20_atomic_functions.h \
    /usr/local/cuda/include/sm_20_atomic_functions.hpp \
    /usr/local/cuda/include/sm_32_atomic_functions.h \
    /usr/local/cuda/include/sm_32_atomic_functions.hpp \
    /usr/local/cuda/include/sm_35_atomic_functions.h \
    /usr/local/cuda/include/sm_60_atomic_functions.h \
    /usr/local/cuda/include/sm_60_atomic_functions.hpp \
    /usr/local/cuda/include/sm_20_intrinsics.h \
    /usr/local/cuda/include/sm_20_intrinsics.hpp \
    /usr/local/cuda/include/sm_30_intrinsics.h \
    /usr/local/cuda/include/sm_30_intrinsics.hpp \
    /usr/local/cuda/include/sm_32_intrinsics.h \
    /usr/local/cuda/include/sm_32_intrinsics.hpp \
    /usr/local/cuda/include/sm_35_intrinsics.h \
    /usr/local/cuda/include/sm_61_intrinsics.h \
    /usr/local/cuda/include/sm_61_intrinsics.hpp \
    /usr/local/cuda/include/crt/sm_70_rt.h \
    /usr/local/cuda/include/crt/sm_70_rt.hpp \
    /usr/local/cuda/include/surface_functions.h \
    /usr/local/cuda/include/texture_fetch_functions.h \
    /usr/local/cuda/include/texture_indirect_functions.h \
    /usr/local/cuda/include/surface_indirect_functions.h \
    /usr/local/cuda/include/device_launch_parameters.h \
    /home/sugimoto/magma/include/magma_v2.h \
    /home/sugimoto/magma/include/magmablas.h \
    /home/sugimoto/magma/include/magma_copy.h \
    /home/sugimoto/magma/include/magma_types.h \
    /usr/local/cuda/include/cuda.h \
    /usr/include/c++/7/bits/stl_relops.h \
    /usr/include/c++/7/bits/stl_pair.h \
    /usr/include/c++/7/initializer_list \
    /usr/local/cuda/include/cusparse_v2.h \
    /usr/local/cuda/include/cusparse.h \
    /home/sugimoto/magma/include/magmablas_z.h \
    /home/sugimoto/magma/include/magmablas_c.h \
    /home/sugimoto/magma/include/magmablas_d.h \
    /home/sugimoto/magma/include/magmablas_s.h \
    /home/sugimoto/magma/include/magmablas_zc.h \
    /home/sugimoto/magma/include/magmablas_ds.h \
    /home/sugimoto/magma/include/magmablas_h.h \
    /home/sugimoto/magma/include/magma_batched.h \
    /home/sugimoto/magma/include/magma_zbatched.h \
    /home/sugimoto/magma/include/magma_cbatched.h \
    /home/sugimoto/magma/include/magma_dbatched.h \
    /home/sugimoto/magma/include/magma_sbatched.h \
    /home/sugimoto/magma/include/magma_hbatched.h \
    /home/sugimoto/magma/include/magma_vbatched.h \
    /home/sugimoto/magma/include/magma_zvbatched.h \
    /home/sugimoto/magma/include/magma_cvbatched.h \
    /home/sugimoto/magma/include/magma_dvbatched.h \
    /home/sugimoto/magma/include/magma_svbatched.h \
    /home/sugimoto/magma/include/magma_bulge.h \
    /home/sugimoto/magma/include/magma_zbulge.h \
    /home/sugimoto/magma/include/magma_cbulge.h \
    /home/sugimoto/magma/include/magma_dbulge.h \
    /home/sugimoto/magma/include/magma_sbulge.h \
    /home/sugimoto/magma/include/magma_z.h \
    /home/sugimoto/magma/include/magma_zgehrd_m.h \
    /home/sugimoto/magma/include/magma_c.h \
    /home/sugimoto/magma/include/magma_cgehrd_m.h \
    /home/sugimoto/magma/include/magma_d.h \
    /home/sugimoto/magma/include/magma_dgehrd_m.h \
    /home/sugimoto/magma/include/magma_s.h \
    /home/sugimoto/magma/include/magma_sgehrd_m.h \
    /home/sugimoto/magma/include/magma_zc.h \
    /home/sugimoto/magma/include/magma_ds.h \
    /home/sugimoto/magma/include/magma_auxiliary.h \
    /home/sugimoto/magma/include/magma_htc.h \
    /home/sugimoto/magma/include/magma_lapack.h \
    /home/sugimoto/magma/include/magma_mangling.h \
    /home/sugimoto/magma/include/magma_mangling_cmake.h \
    /home/sugimoto/magma/include/magma_zlapack.h \
    /home/sugimoto/magma/include/magma_clapack.h \
    /home/sugimoto/magma/include/magma_dlapack.h \
    /home/sugimoto/magma/include/magma_slapack.h \
    /usr/local/cuda/include/curand_kernel.h \
    /usr/local/cuda/include/curand.h \
    /usr/local/cuda/include/curand_discrete.h \
    /usr/local/cuda/include/curand_precalc.h \
    /usr/local/cuda/include/curand_mrg32k3a.h \
    /usr/local/cuda/include/curand_mtgp32_kernel.h \
    /usr/local/cuda/include/curand_mtgp32.h \
    /usr/local/cuda/include/curand_philox4x32_x.h \
    /usr/local/cuda/include/curand_globals.h \
    /usr/local/cuda/include/curand_uniform.h \
    /usr/local/cuda/include/curand_normal.h \
    /usr/local/cuda/include/curand_normal_static.h \
    /usr/local/cuda/include/curand_lognormal.h \
    /usr/local/cuda/include/curand_poisson.h \
    /usr/local/cuda/include/curand_discrete2.h \
    /usr/local/cuda/include/curand_mtgp32_host.h \
    /usr/local/cuda/include/curand_mtgp32dc_p_11213.h
