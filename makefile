all: PeakSeq 

CC = g++
comp_flags = -c -Wall -O3
exec_name = bin/PeakSeq
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/compare_signal_tracks.o \
${LIB_DIR}/enrichment_profile.o \
${LIB_DIR}/input_normalization.o \
${LIB_DIR}/fragment_simulation.o \
${LIB_DIR}/main.o \
${LIB_DIR}/mapped_read_tools.o \
${LIB_DIR}/nomenclature.o \
${LIB_DIR}/annot_region_tools.o \
${LIB_DIR}/peak_region.o \
${LIB_DIR}/peakseq.o \
${LIB_DIR}/rng.o \
${LIB_DIR}/seed_manager.o \
${LIB_DIR}/chip_seq_chr_data.o \
${LIB_DIR}/simulated_background.o \
${LIB_DIR}/poisson_background.o \
${LIB_DIR}/xlog_math.o \
${LIB_DIR}/peakseq_output.o \
${LIB_DIR}/config.o \
${LIB_DIR}/utils.o \
${LIB_DIR}/ansi_string.o \

PeakSeq: ${objs}
	${CC} -O3 -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 
