all: PeakSeq

CC = g++
comp_flags = -c -O3
exec_name = bin/PeakSeq
LIB_DIR = src/lib
FILES_DIR = src/genomics/chip_seq/peak_selector

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/chromosome/chromosome.o \
${FILES_DIR}/compare_signal_tracks.o \
${FILES_DIR}/enrichment_profile.o \
${FILES_DIR}/input_normalization.o \
${FILES_DIR}/main.o \
${LIB_DIR}/database/mapped_read/mapped_read_tools.o \
${LIB_DIR}/database/mapped_read/fragment.o \
${LIB_DIR}/utils/memory/mem_pool.o \
${FILES_DIR}/peak_region.o \
${FILES_DIR}/peakseq.o \
${LIB_DIR}/utils/rng/rng.o \
${LIB_DIR}/utils/rng/seed_manager.o \
${FILES_DIR}/chip_seq_chr_data.o \
${FILES_DIR}/simulated_background.o \
${FILES_DIR}/poisson_background.o \
${LIB_DIR}/utils/xmath/log/xlog_math.o \
${LIB_DIR}/utils/xmath/linear/linear_math.o \
${FILES_DIR}/peakseq_output.o \
${LIB_DIR}/utils/ansi_cli/config.o \
${LIB_DIR}/utils/file/utils.o \
${LIB_DIR}/utils/ansi_string/ansi_string.o \

PeakSeq: ${objs}
	${CC} -O3 -o ${exec_name} ${objs}

clean:
	find . -name *.o | xargs -Ifiles rm -f files
