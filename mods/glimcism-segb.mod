*DECLARE makefileu_comp_in
*I makefileu_comp_in.5
ICE_DIR=/storage/basic/baobab/seg/FAMOUS
GLIMMER=$(ICE_DIR)/UNIC9is/trunk/cism-parallel
NETCDF=$(ICE_DIR)/NCDIRF2

VPATH=$(GLIMMER)/include
*DECLARE makefileu_link_in
*I makefileu_link_in.5
GCOM3=$(UMDIR)/gcom3.0

ICE_DIR=/storage/basic/baobab/seg/FAMOUS
HDF_DIR=$(ICE_DIR)/NCDIRF2
CDF_DIR=$(ICE_DIR)/NCDIRF2
Z_DIR=$(ICE_DIR)/NCDIRF2
GLIM_DIR=$(HOME)/Data/plane-17/FAMOUS/UNIC9is/trunk/cism-parallel
GLIM_DIR=$(ICE_DIR)/UNIC9is/trunk/cism-parallel
CDF_LIBS=-L$(CDF_DIR)/lib -L$(HDF_DIR)/lib -L$(Z_DIR)/lib \
         -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
GLIM_LIBS=-L$(GLIM_DIR)/lib -lglint -lglide -lglimmer \
          -lglimmer-solve -lglimmer-IO

*DECLARE qsmncompile
*D PXCLLMOD.50
${MAKEFILEUCOMPIN:-$TEMP/modscr_$RUNID/makefileu_comp_in}
*D PXCLLMOD.52
${MAKEFILEULINKIN:-$TEMP/modscr_$RUNID/makefileu_link_in}
