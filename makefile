MRT_INSTALL = /opt/MRT

CC = gcc
CFLAGS = -g
LDFLAGS = $(GEOLIB) -lm -lz -s

CP = cp
MV = mv
RM = rm -f
INCS = -I$(MRT_INSTALL)/include -I$(MRT_INSTALL)/gctp -I$(MRT_INSTALL)/geolib


SRC     = HLS_VIIRS_shape_fusion.c

OBJ = $(SRC:.c=.o)


GEOLIB = $(MRT_INSTALL)/geolib/libgeolib.a $(MRT_INSTALL)/gctp/libgctp.a


PLM_LSPD_07_HLS_VIIRS_4band: $(SRC)
	$(CC) $(SRC) $(CFLAGS) $(INCS) $(LDFLAGS) -o HLS_VIIRS_shape_fusion.exe




