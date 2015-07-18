SET( PGVL_SRCS
   Image.cpp
   ImageProcessing.cpp
   ppm.cpp
)

ADD_LIBRARY( pgvl
   SHARED
   ${PGVL_SRCS}
)
